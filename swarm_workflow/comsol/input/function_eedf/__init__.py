"""Canonical moment-conserving Function-EEDF representations for COMSOL."""

from swarm_workflow.tables.contracts import (
    COLLISION_RATE_KERNEL_COLUMNS,
    COLLISION_RATE_KERNEL_TABLE,
)

from .contracts import (
    COMSOL_EEDF_COLUMNS,
    SOURCE_MEAN_RELATIVE_ERROR_LIMIT,
    C1FunctionEedf,
    CollisionRateKernel,
    ComsolEedfImportContract,
    ComsolFunctionEedfGrid,
    FunctionEedfError,
)
from .c1 import (
    MEAN_AXIS_RATE_IMPORTANCE_FRACTION,
    RATE_IMPORTANCE_FRACTION,
    RATE_SCALED_ERROR_TOLERANCE,
    SHAPE_TOTAL_VARIATION_TOLERANCE,
    build_c1_function_eedf,
    project_c1_function_eedf_to_comsol_grid,
)
from .io import read_c1_function_eedf, read_comsol_function_eedf_grid
from .kernels import collision_rate_coefficient, read_collision_rate_kernels
from .moments import (
    evaluate_c1_function_eedf,
    evaluate_comsol_function_eedf_grid,
    pchip_weighted_moments,
    piecewise_linear_weighted_moments,
    scaled_shape_moments,
)

__all__ = (
    "COMSOL_EEDF_COLUMNS",
    "MEAN_AXIS_RATE_IMPORTANCE_FRACTION",
    "SOURCE_MEAN_RELATIVE_ERROR_LIMIT",
    "RATE_IMPORTANCE_FRACTION",
    "RATE_SCALED_ERROR_TOLERANCE",
    "SHAPE_TOTAL_VARIATION_TOLERANCE",
    "C1FunctionEedf",
    "COLLISION_RATE_KERNEL_COLUMNS",
    "COLLISION_RATE_KERNEL_TABLE",
    "CollisionRateKernel",
    "ComsolEedfImportContract",
    "ComsolFunctionEedfGrid",
    "FunctionEedfError",
    "build_c1_function_eedf",
    "collision_rate_coefficient",
    "evaluate_c1_function_eedf",
    "evaluate_comsol_function_eedf_grid",
    "pchip_weighted_moments",
    "piecewise_linear_weighted_moments",
    "project_c1_function_eedf_to_comsol_grid",
    "read_c1_function_eedf",
    "read_collision_rate_kernels",
    "read_comsol_function_eedf_grid",
    "scaled_shape_moments",
)
