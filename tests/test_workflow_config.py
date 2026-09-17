from __future__ import annotations

import pytest

import swarm_workflow
import swarm_workflow.campaign.sweep as sweep_module
from swarm_workflow.campaign.config import (
    DeterministicExecutionConfig,
    MeanEnergySupportConfig,
    MixtureSpec,
    WorkflowConfig,
    load_workflow,
)


def test_package_exports_workflow_config_from_its_owner() -> None:
    assert swarm_workflow.DeterministicExecutionConfig is DeterministicExecutionConfig
    assert swarm_workflow.MeanEnergySupportConfig is MeanEnergySupportConfig
    assert swarm_workflow.MixtureSpec is MixtureSpec
    assert swarm_workflow.WorkflowConfig is WorkflowConfig
    assert swarm_workflow.load_workflow is load_workflow


@pytest.mark.parametrize(
    "name",
    (
        "DeterministicExecutionConfig",
        "MeanEnergySupportConfig",
        "MixtureSpec",
        "WorkflowConfig",
        "load_workflow",
        "_combined_hash",
        "_cross_section_file_hashes",
        "_sha256_file",
    ),
)
def test_sweep_does_not_reexport_workflow_config_api(name: str) -> None:
    assert not hasattr(sweep_module, name)
