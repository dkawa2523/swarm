"""Public API for provenance-bound COMSOL Function-EEDF audits."""

from .analysis import analyze_comsol_eedf_audit
from .contracts import ComsolEedfAuditError, ComsolEedfAuditPlan
from .io import read_native_eedf_grid
from .java import extract_comsol_eedf_audit_log, render_comsol_eedf_audit_java
from .planning import prepare_comsol_eedf_audit

__all__ = [
    "ComsolEedfAuditError",
    "ComsolEedfAuditPlan",
    "analyze_comsol_eedf_audit",
    "extract_comsol_eedf_audit_log",
    "prepare_comsol_eedf_audit",
    "read_native_eedf_grid",
    "render_comsol_eedf_audit_java",
]
