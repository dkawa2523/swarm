param(
    [Parameter(Mandatory = $true)]
    [int]$DatabaseId,

    [Parameter(Mandatory = $true)]
    [int]$TargetSpeciesId,

    [Parameter(Mandatory = $true)]
    [string]$OutputPath,

    [string]$Mirror = "https://nl.lxcat.net",

    [string[]]$ExcludeProcessType = @(),

    [string[]]$ExcludeProcessLabelRegex = @(),

    [switch]$Force
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$resolvedOutput = [System.IO.Path]::GetFullPath($OutputPath)
if ((Test-Path -LiteralPath $resolvedOutput) -and -not $Force) {
    throw "Output already exists: $resolvedOutput (pass -Force to replace it)"
}

$parent = Split-Path -Parent $resolvedOutput
if (-not (Test-Path -LiteralPath $parent)) {
    New-Item -ItemType Directory -Path $parent | Out-Null
}

$session = [Microsoft.PowerShell.Commands.WebRequestSession]::new()
$selectionUri = "$Mirror/data/preselect.php?t=cs&d=$DatabaseId&a=1&b=$TargetSpeciesId"
$selection = Invoke-WebRequest -Uri $selectionUri -WebSession $session -UseBasicParsing
$processRows = [regex]::Matches(
    $selection.Content,
    '(?s)<li><input[^>]*name="proc\[\]" value="(\d+)"[^>]*>.*?<label[^>]*>(.*?)</label></li>'
) | ForEach-Object {
    $label = [regex]::Replace($_.Groups[2].Value, '<[^>]+>', '')
    [pscustomobject]@{
        id = $_.Groups[1].Value
        label = [System.Net.WebUtility]::HtmlDecode($label).Trim()
    }
}

$excludedTypes = @($ExcludeProcessType | ForEach-Object { $_.Trim().ToLowerInvariant() })
$excludedLabelPatterns = @($ExcludeProcessLabelRegex | Where-Object { $_.Length -gt 0 })
$selectedRows = @($processRows | Where-Object {
    $firstWord = ($_.label -split '\s+', 2)[0].ToLowerInvariant()
    $labelExcluded = $false
    foreach ($pattern in $excludedLabelPatterns) {
        if ($_.label -match $pattern) {
            $labelExcluded = $true
            break
        }
    }
    ($firstWord -notin $excludedTypes) -and -not $labelExcluded
})
$excludedRows = @($processRows | Where-Object { $_.id -notin $selectedRows.id })
$processIds = @($selectedRows | ForEach-Object { $_.id })

if ($processIds.Count -eq 0) {
    throw "LXCat returned no selectable processes for $selectionUri"
}

$postBody = ($processIds | ForEach-Object {
    "proc%5B%5D=$([uri]::EscapeDataString($_))"
}) -join "&"

function Submit-Selection {
    Invoke-WebRequest `
        -Uri "$Mirror/data/set_processes.php" `
        -WebSession $session `
        -Method Post `
        -ContentType "application/x-www-form-urlencoded" `
        -Body $postBody `
        -UseBasicParsing
}

$response = Submit-Selection
$acceptMatch = [regex]::Match(
    $response.Content,
    'class="accept"><a href="([^"]+)"'
)
if ($acceptMatch.Success) {
    $acceptUri = [uri]::new([uri]$Mirror, $acceptMatch.Groups[1].Value)
    Invoke-WebRequest -Uri $acceptUri -WebSession $session -UseBasicParsing | Out-Null
    $response = Submit-Selection
}

$jobUri = $response.BaseResponse.RequestMessage.RequestUri.AbsoluteUri
$outputPage = $null
for ($attempt = 1; $attempt -le 30; $attempt++) {
    Start-Sleep -Milliseconds 500
    $candidate = Invoke-WebRequest `
        -Uri $jobUri `
        -WebSession $session `
        -UseBasicParsing `
        -AllowInsecureRedirect
    if ($candidate.Content -match '<title>output</title>') {
        $outputPage = $candidate
        break
    }
}

if ($null -eq $outputPage) {
    throw "LXCat did not finish preparing the selected data within 15 seconds"
}

$cacheUri = $outputPage.BaseResponse.RequestMessage.RequestUri
$downloadUri = [uri]::new($cacheUri, "Cross%20section.txt")
Invoke-WebRequest `
    -Uri $downloadUri `
    -WebSession $session `
    -UseBasicParsing `
    -AllowInsecureRedirect `
    -OutFile $resolvedOutput

$header = Get-Content -LiteralPath $resolvedOutput -TotalCount 8
if (($header -join "`n") -notmatch '^LXCat, www\.lxcat\.net') {
    throw "Downloaded file is not an LXCat legacy cross-section file: $resolvedOutput"
}

$hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $resolvedOutput).Hash.ToLowerInvariant()
[pscustomobject]@{
    selection_url = $selectionUri
    downloaded_url = $downloadUri.AbsoluteUri
    database_id = $DatabaseId
    target_species_id = $TargetSpeciesId
    excluded_process_types = $excludedTypes
    excluded_process_label_regex = $excludedLabelPatterns
    excluded_process_labels = @($excludedRows | ForEach-Object { $_.label })
    available_process_count = $processRows.Count
    process_count = $processIds.Count
    retrieved_at_utc = [DateTime]::UtcNow.ToString("o")
    output_path = $resolvedOutput
    sha256 = $hash
} | ConvertTo-Json
