"""Shared Java literal helpers for COMSOL source generators."""

from __future__ import annotations

import re
from collections.abc import Sequence


_PUBLIC_CLASS = re.compile(r"\bpublic\s+class\s+([A-Za-z_$][A-Za-z0-9_$]*)\b")


def java_path(value: object) -> str:
    """Return a path string that is stable inside generated Java source."""

    return str(value).replace("\\", "/")


def java_string(value: object) -> str:
    """Return ``value`` as an escaped Java string literal."""

    text = str(value)
    escaped = (
        text.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\n", "\\n")
        .replace("\r", "\\r")
    )
    return f'"{escaped}"'


def compose_java_mains(
    class_name: str,
    stages: Sequence[tuple[str, str]],
) -> str:
    """Compose generated Java programs into one ordered COMSOL entry point.

    The input sources are repository-generated, package-free Java programs
    with one public class and a conventional ``main`` method.  Package-private
    copies let COMSOL compile the public entry point with all symbols resolved;
    callers must also compile and stage the original public sources beside it
    because ``comsolcompile`` emits only the requested public class.  Keeping
    each stage as its own generated source leaves model-specific generators
    independently testable while avoiding a new COMSOL JVM for each read-only
    verification stage.
    """

    if not re.fullmatch(r"[A-Za-z_$][A-Za-z0-9_$]*", class_name):
        raise ValueError(f"invalid Java class name: {class_name!r}")
    if not stages:
        raise ValueError("at least one Java stage is required")

    imports: list[str] = []
    bodies: list[str] = []
    seen_classes: set[str] = set()
    for expected_class, source in stages:
        if expected_class in seen_classes:
            raise ValueError(f"duplicate Java stage class: {expected_class}")
        seen_classes.add(expected_class)
        if re.search(r"(?m)^\s*package\s+", source):
            raise ValueError("generated COMSOL Java stages must not declare a package")
        matches = _PUBLIC_CLASS.findall(source)
        if matches != [expected_class]:
            raise ValueError(
                "generated Java stage must contain exactly its declared public class: "
                f"expected {expected_class!r}, found {matches!r}"
            )
        if "public static void main(" not in source:
            raise ValueError(
                f"generated Java stage has no public main: {expected_class}"
            )

        body_lines: list[str] = []
        for line in source.splitlines():
            stripped = line.strip()
            if stripped.startswith("import "):
                if stripped not in imports:
                    imports.append(stripped)
                continue
            body_lines.append(line)
        body = "\n".join(body_lines).strip()
        body = re.sub(
            rf"\bpublic\s+class\s+{re.escape(expected_class)}\b",
            f"final class {expected_class}",
            body,
            count=1,
        )
        bodies.append(body)

    calls = "\n".join(f"    {stage}.main(args);" for stage, _ in stages)
    wrapper = (
        f"public class {class_name} {{\n"
        "  public static void main(String[] args) throws Exception {\n"
        f"{calls}\n"
        "  }\n"
        "}"
    )
    return "\n\n".join(("\n".join(imports), *bodies, wrapper, ""))


__all__ = ["compose_java_mains", "java_path", "java_string"]
