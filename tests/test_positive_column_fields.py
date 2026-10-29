"""Configured saved-solution selection produces actual, identifiable CSV rows."""

import csv
import io
from dataclasses import replace
import os
from pathlib import Path
from swarm_workflow.comsol.models.positive_column.results import read_mapped_results
import shutil
import subprocess

import pytest
import yaml

from swarm_workflow.comsol.models.positive_column.config import (
    ComsolMappingError,
    ResultField,
    load_comsol_mapping,
)
from swarm_workflow.comsol.models.positive_column.java import (
    generate_field_export_java_source,
)
from test_comsol_positive_column import _write_positive_column_repo
from support.positive_column import result_csv
from swarm_workflow.comsol.models.positive_column.audit import (
    validate_positive_column_result,
)
from swarm_workflow.comsol.models.positive_column.workflow import (
    prepare_positive_column_run,
    PositiveColumnWorkflowError,
)


@pytest.mark.parametrize(
    "mutation,reason",
    [
        ({"probes": []}, "migrate"),
        ({"point_m": [0.1]}, "unsupported results"),
        ({"solutions": "automatic"}, "all or last"),
        (
            {
                "fields": [
                    {
                        "name": "ne",
                        "expression": "plas.ne",
                        "unit": "1/m^3",
                        "time_s": 0.1,
                    }
                ]
            },
            "unsupported results.fields",
        ),
        (
            {"fields": [{"name": "ne", "expression": "plas.ne", "unit": "1/m^3"}] * 2},
            "duplicate",
        ),
        ({"fields": [{"name": "time_s", "expression": "t", "unit": "s"}]}, "reserved"),
        (
            {
                "fields": [
                    {
                        "name": "mean_electron_energy",
                        "expression": "plas.te",
                        "unit": "K",
                    }
                ]
            },
            "unit.*unsupported",
        ),
    ],
)
def test_invalid_result_requests_fail_instead_of_changing_meaning(
    tmp_path, mutation, reason
):
    path, _ = _write_positive_column_repo(tmp_path)
    raw = yaml.safe_load(path.read_text())
    raw["results"].update(mutation)
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ComsolMappingError, match=reason):
        load_comsol_mapping(path)


@pytest.mark.parametrize(
    "fault,reason",
    [
        ("dataset", "dataset differs"),
        ("duplicate", "duplicate solution/point"),
        ("multiple", "one saved solution"),
        ("time", "mixes parameter/time"),
    ],
)
def test_result_acceptance_requires_the_selected_solution_identity(
    tmp_path, fault, reason
):
    path, bundle = _write_positive_column_repo(tmp_path)
    plan = prepare_positive_column_run(path, bundle_path=bundle)
    text = result_csv(
        "E_over_N,mean_electron_energy,reduced_mobility,reduced_diffusion_L,excitation_townsend,ionization_townsend\n"
        "4e-20,4,10.6666666666667,20.6666666666667,1.333333333333e-21,1.333333333333e-22\n"
    )
    reader = csv.DictReader(io.StringIO(text))
    columns = reader.fieldnames
    rows = list(reader)
    if fault == "dataset":
        rows[0]["dataset"] = "unselected_dataset"
    else:
        rows.append(dict(rows[0]))
        if fault == "multiple":
            rows[1]["solution_index"] = "2"
        if fault == "time":
            rows[1]["point_index"] = "1"
            for index, row in enumerate(rows):
                row.update(
                    parameter_name="t",
                    parameter_unit="s",
                    parameter_value=str(index),
                    time_s=str(index),
                )
    with plan.result_csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, columns)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(PositiveColumnWorkflowError, match=reason):
        validate_positive_column_result(
            plan, read_mapped_results(plan.result_csv_path, plan.mapping.results)
        )


@pytest.fixture(scope="module")
def java_api(tmp_path_factory):
    directory = tmp_path_factory.mktemp("field-java-api")
    candidate = (
        Path(os.environ.get("JAVA_HOME", ""))
        / "bin"
        / ("java.exe" if os.name == "nt" else "java")
    )
    java = str(candidate) if candidate.is_file() else shutil.which("java")
    if java is None:
        pytest.skip("generated Java behavior requires a JDK compiler runtime")
    sources = {
        "com/comsol/model/SolutionInfo.java": """package com.comsol.model;
public class SolutionInfo {
  public int[] getOuterSolnum() { return Boolean.getBoolean("outer") ? new int[]{1,2} : new int[]{1}; }
  public int[] getSolnum(int outer, boolean strict) { return new int[]{1,2}; }
  public String getUnit(String name) { return "s"; }
}""",
        "com/comsol/model/util/ModelUtil.java": """package com.comsol.model.util;
import com.comsol.model.Model;
public class ModelUtil {
  public static Model loadCopy(String name, String file) { return new Model(); }
}""",
        "com/comsol/model/Model.java": """package com.comsol.model;
public class Model {
  private final Result result = new Result();
  public Result result() { return result; }
  public Solution sol(String name) { if (!name.equals("sol1")) throw new RuntimeException("wrong solution"); return new Solution(); }
  public static class Solution {
    public SolutionInfo getSolutioninfo() { return new SolutionInfo(); }
    public String[] getPNames() { return new String[]{"t"}; }
    public double[] getPVals() { return new double[]{0.001,0.002}; }
  }
  public static class Result {
    private final Numerical numerical = new Numerical();
    public Dataset dataset(String name) { if (!name.equals("dset1")) throw new RuntimeException("wrong dataset"); return new Dataset(); }
    public NumericalList numerical() { return new NumericalList(); }
    public Numerical numerical(String name) { return numerical; }
  }
  public static class Dataset { public String getString(String name) { return "sol1"; } }
  public static class NumericalList { public void remove(String name) {} public void create(String name,String type) {} }
  public static class Numerical {
    private String dataset;
    private String[] expressions;
    private String[] units;
    public void set(String name, String value) { if (name.equals("data")) dataset=value; }
    public void set(String name, String[] values) { if (name.equals("expr")) expressions=values; else if (name.equals("unit")) units=values; }
    public String[] getStringArray(String name) { String[] result=units.clone(); if (Boolean.getBoolean("unitMismatch")) result[0]="cm"; return result; }
    public boolean isComplex() { return Boolean.getBoolean("complex"); }
    public int getNData() { return 3; }
    public double[][][] getData() {
      if (!"dset1".equals(dataset)) throw new RuntimeException("dataset was not selected");
      double[][][] result=new double[expressions.length][2][3];
      for (int e=0;e<expressions.length;e++) for (int s=0;s<2;s++) for (int p=0;p<3;p++) {
        String name=expressions[e];
        if (name.equals("x")) result[e][s][p]=0.01*p;
        else if (name.equals("plas.ne")) result[e][s][p]=1000*(s+1)+p;
        else if (name.equals("t")) result[e][s][p]=0.001*(s+1);
        else throw new RuntimeException("unknown expression");
      }
      if (Boolean.getBoolean("nonfinite")) result[0][0][0]=Double.NaN;
      return result;
    }
  }
}""",
    }
    for name, source in sources.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(source)
    classes = directory / "classes"
    command = [
        java,
        "-m",
        "jdk.compiler/com.sun.tools.javac.Main",
        "-d",
        str(classes),
        *[str(directory / name) for name in sources],
    ]
    result = subprocess.run(
        command, capture_output=True, text=True, encoding="utf-8", errors="replace"
    )
    assert result.returncode == 0, result.stderr
    return java, classes


@pytest.mark.parametrize(
    "selection,expected_solutions", [("all", [1, 2]), ("last", [2])]
)
def test_generated_export_writes_selected_solutions_and_times(
    tmp_path, java_api, selection, expected_solutions
):
    java, classes = java_api
    path, _ = _write_positive_column_repo(tmp_path)
    mapping = load_comsol_mapping(path)
    mapping = replace(
        mapping,
        results=replace(
            mapping.results,
            solutions=selection,
            fields=(
                ResultField("x", "x", "m"),
                ResultField("ne", "plas.ne", "1/m^3"),
            ),
        ),
    )
    output = tmp_path / "fields.csv"
    source = tmp_path / "FieldExport.java"
    source.write_text(
        generate_field_export_java_source(
            mapping, class_name="FieldExport", output_csv=output
        )
    )
    compile_result = subprocess.run(
        [
            java,
            "-m",
            "jdk.compiler/com.sun.tools.javac.Main",
            "-cp",
            str(classes),
            "-d",
            str(tmp_path),
            str(source),
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert compile_result.returncode == 0, compile_result.stderr
    classpath = os.pathsep.join((str(tmp_path), str(classes)))
    executed = subprocess.run(
        [java, "-cp", classpath, "FieldExport"], capture_output=True, text=True
    )
    assert executed.returncode == 0, executed.stderr
    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3 * len(expected_solutions)
    assert sorted({int(row["solution_index"]) for row in rows}) == expected_solutions
    for row in rows:
        solution, point = int(row["solution_index"]), int(row["point_index"])
        assert row["dataset"] == "dset1"
        assert float(row["time_s"]) == pytest.approx(0.001 * solution)
        assert row["parameter_name"] == "t" and row["parameter_unit"] == "s"
        assert float(row["parameter_value"]) == pytest.approx(float(row["time_s"]))
        assert float(row["x"]) == pytest.approx(0.01 * point)
        assert float(row["ne"]) == 1000 * solution + point
    original = output.read_bytes()
    for option, reason in (
        ("-Douter=true", "outer parameter sweeps"),
        ("-Dnonfinite=true", "nonfinite"),
        ("-DunitMismatch=true", "evaluated unit differs"),
        ("-Dcomplex=true", "refusing to discard imaginary"),
    ):
        failed = subprocess.run(
            [java, option, "-cp", classpath, "FieldExport"],
            capture_output=True,
            text=True,
        )
        assert failed.returncode != 0 and reason in failed.stderr
        assert output.read_bytes() == original
