"""Public model-independent COMSOL execution surface."""

from .commands import build_java_batch_command, build_java_compile_command
from .executable import resolve_comsol_executable
from .runner import execute_comsol_java_source, execute_generated_comsol_java
from .types import (
    ComsolAdapterError,
    ComsolExecutionSummary,
    ComsolRunContext,
    ResolvedComsolExecutable,
)

__all__ = (
    "ComsolAdapterError",
    "ComsolExecutionSummary",
    "ComsolRunContext",
    "ResolvedComsolExecutable",
    "build_java_batch_command",
    "build_java_compile_command",
    "execute_comsol_java_source",
    "execute_generated_comsol_java",
    "resolve_comsol_executable",
)
