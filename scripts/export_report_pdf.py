from __future__ import annotations

import html
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT_MD = ROOT / "report.md"
OUT_PDF = ROOT / "report.pdf"
PDF_ASSET_DIR = ROOT / "outputs" / "ar_n2_ratio_temp_grid" / "plots" / "pdf_assets"


def register_fonts() -> tuple[str, str]:
    font = "HeiseiKakuGo-W5"
    mono = "Courier"
    pdfmetrics.registerFont(UnicodeCIDFont(font))
    return font, mono


def make_styles(font: str, mono: str):
    base = getSampleStyleSheet()
    styles = {
        "title": ParagraphStyle(
            "JapaneseTitle",
            parent=base["Title"],
            fontName=font,
            fontSize=20,
            leading=26,
            alignment=TA_CENTER,
            spaceAfter=14,
        ),
        "h1": ParagraphStyle(
            "JapaneseH1",
            parent=base["Heading1"],
            fontName=font,
            fontSize=16,
            leading=21,
            spaceBefore=10,
            spaceAfter=8,
        ),
        "h2": ParagraphStyle(
            "JapaneseH2",
            parent=base["Heading2"],
            fontName=font,
            fontSize=13,
            leading=18,
            spaceBefore=8,
            spaceAfter=6,
        ),
        "h3": ParagraphStyle(
            "JapaneseH3",
            parent=base["Heading3"],
            fontName=font,
            fontSize=11.5,
            leading=16,
            spaceBefore=6,
            spaceAfter=4,
        ),
        "body": ParagraphStyle(
            "JapaneseBody",
            parent=base["BodyText"],
            fontName=font,
            fontSize=9.2,
            leading=14,
            alignment=TA_LEFT,
            spaceAfter=5,
        ),
        "caption": ParagraphStyle(
            "JapaneseCaption",
            parent=base["BodyText"],
            fontName=font,
            fontSize=8.3,
            leading=11,
            alignment=TA_CENTER,
            textColor=colors.HexColor("#444444"),
            spaceBefore=2,
            spaceAfter=8,
        ),
        "code": ParagraphStyle(
            "Code",
            parent=base["Code"],
            fontName=mono,
            fontSize=7.2,
            leading=9,
            leftIndent=4,
            rightIndent=4,
            borderPadding=6,
            backColor=colors.HexColor("#f4f4f4"),
            borderColor=colors.HexColor("#dddddd"),
            borderWidth=0.5,
            spaceBefore=3,
            spaceAfter=7,
        ),
        "table": ParagraphStyle(
            "JapaneseTable",
            parent=base["BodyText"],
            fontName=font,
            fontSize=7.2,
            leading=9,
        ),
        "table_head": ParagraphStyle(
            "JapaneseTableHead",
            parent=base["BodyText"],
            fontName=font,
            fontSize=7.4,
            leading=9.2,
            textColor=colors.white,
        ),
    }
    return styles


def clean_inline(text: str, mono: str = "Courier") -> str:
    escaped = html.escape(text.strip())

    def repl(match: re.Match[str]) -> str:
        return f'<font name="{mono}">{match.group(1)}</font>'

    escaped = re.sub(r"`([^`]+)`", repl, escaped)
    escaped = escaped.replace(" x ", " &times; ")
    return escaped


def flush_paragraph(buf: list[str], story: list, styles) -> None:
    if not buf:
        return
    text = " ".join(line.strip() for line in buf).strip()
    if text:
        story.append(Paragraph(clean_inline(text), styles["body"]))
    buf.clear()


def split_table_row(line: str) -> list[str]:
    stripped = line.strip().strip("|")
    return [cell.strip() for cell in stripped.split("|")]


def is_table_separator(line: str) -> bool:
    return bool(re.match(r"^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)+\|?\s*$", line))


def add_table(rows: list[list[str]], story: list, styles, available_width: float) -> None:
    if len(rows) < 2:
        return
    header = rows[0]
    body = [row for row in rows[1:] if not all(re.fullmatch(r":?-{3,}:?", c) for c in row)]
    ncols = max(len(header), *(len(r) for r in body)) if body else len(header)
    normalized = []
    for row in [header] + body:
        padded = row + [""] * (ncols - len(row))
        normalized.append(padded[:ncols])

    data = []
    for ridx, row in enumerate(normalized):
        style = styles["table_head"] if ridx == 0 else styles["table"]
        data.append([Paragraph(clean_inline(cell), style) for cell in row])

    col_width = available_width / ncols
    table = Table(data, colWidths=[col_width] * ncols, repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#334155")),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#cbd5e1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f8fafc")]),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 0.22 * cm))


def add_image(markdown_line: str, story: list, styles, available_width: float) -> None:
    match = re.match(r"!\[([^\]]*)\]\(([^)]+)\)", markdown_line.strip())
    if not match:
        return
    alt, rel_path = match.groups()
    img_path = (ROOT / rel_path).resolve()
    if not img_path.exists():
        story.append(Paragraph(f"Missing image: {html.escape(rel_path)}", styles["body"]))
        return

    with PILImage.open(img_path) as im:
        width_px, height_px = im.size
    draw_width = min(available_width, 16.5 * cm)
    draw_height = draw_width * height_px / width_px
    if draw_height > 17.5 * cm:
        draw_height = 17.5 * cm
        draw_width = draw_height * width_px / height_px

    block = [
        Image(str(img_path), width=draw_width, height=draw_height),
        Paragraph(clean_inline(alt), styles["caption"]),
    ]
    story.append(KeepTogether(block))


def render_math_image(math_text: str, index: int) -> Path:
    PDF_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    eq = " ".join(line.strip() for line in math_text.splitlines() if line.strip())
    eq = eq.replace(r"\boldsymbol", r"\mathbf")

    width = min(max(5.8, len(eq) * 0.07), 10.8)
    fig = plt.figure(figsize=(width, 0.85))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.5, f"${eq}$", ha="center", va="center", fontsize=15)
    out = PDF_ASSET_DIR / f"equation_{index:02d}.png"
    fig.savefig(out, dpi=220, bbox_inches="tight", pad_inches=0.08, transparent=False)
    plt.close(fig)
    return out


def add_math_block(math_text: str, story: list, styles, available_width: float, index: int) -> None:
    img_path = render_math_image(math_text, index)
    with PILImage.open(img_path) as im:
        width_px, height_px = im.size
    draw_width = min(available_width * 0.9, width_px * 0.72)
    draw_height = draw_width * height_px / width_px
    story.append(KeepTogether([Image(str(img_path), width=draw_width, height=draw_height)]))
    story.append(Spacer(1, 0.12 * cm))


def add_code_block(lang: str, code: list[str], story: list, styles, available_width: float, math_index: int) -> int:
    label = None
    if lang == "mermaid":
        label = "Mermaid workflow definition"
    elif lang == "math":
        label = "TeX equation"
    if label:
        story.append(Paragraph(label, styles["h3"]))
    text = "\n".join(code).strip()
    story.append(Preformatted(text, styles["code"], maxLineLength=92))
    return math_index


def build_story(markdown: str, styles, available_width: float) -> list:
    lines = markdown.splitlines()
    story: list = []
    para: list[str] = []
    table_rows: list[list[str]] = []
    in_code = False
    in_tex = False
    code_lang = ""
    code_lines: list[str] = []
    tex_lines: list[str] = []
    first_heading = True
    math_index = 1
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if in_code:
            if stripped.startswith("```"):
                math_index = add_code_block(
                    code_lang, code_lines, story, styles, available_width, math_index
                )
                in_code = False
                code_lang = ""
                code_lines = []
            else:
                code_lines.append(line)
            i += 1
            continue

        if in_tex:
            if stripped == "$$":
                add_math_block(
                    "\n".join(tex_lines).strip(),
                    story,
                    styles,
                    available_width,
                    math_index,
                )
                math_index += 1
                in_tex = False
                tex_lines = []
            else:
                tex_lines.append(line)
            i += 1
            continue

        if stripped == "$$":
            flush_paragraph(para, story, styles)
            if table_rows:
                add_table(table_rows, story, styles, available_width)
                table_rows = []
            in_tex = True
            tex_lines = []
            i += 1
            continue

        if stripped.startswith("```"):
            flush_paragraph(para, story, styles)
            if table_rows:
                add_table(table_rows, story, styles, available_width)
                table_rows = []
            code_lang = stripped[3:].strip().lower()
            in_code = True
            i += 1
            continue

        if stripped.startswith("|") and "|" in stripped[1:]:
            flush_paragraph(para, story, styles)
            if not is_table_separator(stripped):
                table_rows.append(split_table_row(stripped))
            i += 1
            continue
        elif table_rows:
            add_table(table_rows, story, styles, available_width)
            table_rows = []

        if not stripped:
            flush_paragraph(para, story, styles)
            i += 1
            continue

        if stripped.startswith("!["):
            flush_paragraph(para, story, styles)
            add_image(stripped, story, styles, available_width)
            i += 1
            continue

        if stripped.startswith("#"):
            flush_paragraph(para, story, styles)
            level = len(stripped) - len(stripped.lstrip("#"))
            title = stripped[level:].strip()
            if level == 1:
                if first_heading:
                    story.append(Paragraph(clean_inline(title), styles["title"]))
                    first_heading = False
                else:
                    story.append(PageBreak())
                    story.append(Paragraph(clean_inline(title), styles["h1"]))
            elif level == 2:
                story.append(Paragraph(clean_inline(title), styles["h1"]))
            elif level == 3:
                story.append(Paragraph(clean_inline(title), styles["h2"]))
            else:
                story.append(Paragraph(clean_inline(title), styles["h3"]))
            i += 1
            continue

        if stripped.startswith("- "):
            flush_paragraph(para, story, styles)
            bullet = stripped[2:].strip()
            story.append(Paragraph(f"&bull; {clean_inline(bullet)}", styles["body"]))
            i += 1
            continue

        para.append(line)
        i += 1

    flush_paragraph(para, story, styles)
    if table_rows:
        add_table(table_rows, story, styles, available_width)
    return story


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("HeiseiKakuGo-W5", 7.5)
    canvas.setFillColor(colors.HexColor("#64748b"))
    canvas.drawString(doc.leftMargin, 0.8 * cm, "Swarm solver project notes")
    canvas.drawRightString(A4[0] - doc.rightMargin, 0.8 * cm, f"{doc.page}")
    canvas.restoreState()


def main() -> None:
    font, mono = register_fonts()
    styles = make_styles(font, mono)
    doc = SimpleDocTemplate(
        str(OUT_PDF),
        pagesize=A4,
        rightMargin=1.35 * cm,
        leftMargin=1.35 * cm,
        topMargin=1.3 * cm,
        bottomMargin=1.25 * cm,
        title="Swarm solver project notes",
        author="swarm",
    )
    markdown = REPORT_MD.read_text(encoding="utf-8")
    story = build_story(markdown, styles, doc.width)
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
