"""Sample qualified DC Ar coefficients for a model at 300 K and 2.66644 Pa."""

import argparse
import json
from pathlib import Path

from electron_swarm.core.constants import BOLTZMANN_J_K
from swarm_workflow.tables import Reaction, load_coefficients


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", type=Path, help="one build-tables mixture directory")
    parser.add_argument("--e-over-n-Td", type=float, required=True)
    args = parser.parse_args()

    temperature_K = 300.0
    coefficients = load_coefficients(
        args.tables,
        gas_temperature_K=temperature_K,
        gas_number_density_m3=2.66644 / (BOLTZMANN_J_K * temperature_K),
        mole_fractions={"Ar": 1.0},
    )
    ionization = coefficients.series(
        "rate_coefficient_m3_s",
        reaction=Reaction("Ar", "e+Ar=>2e+Ar+", "ionization", 15.8),
    )
    mobility = coefficients.series("reduced_mobility_m2_V_s_m3")
    print(
        json.dumps(
            {
                "E_over_N_Td": args.e_over_n_Td,
                "source": coefficients.source,
                "ionization": {
                    "value": ionization.at(args.e_over_n_Td),
                    "unit": ionization.unit,
                },
                "reduced_mobility": {
                    "value": mobility.at(args.e_over_n_Td),
                    "unit": mobility.unit,
                },
                "source_manifest_sha256": coefficients.manifest_sha256,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
