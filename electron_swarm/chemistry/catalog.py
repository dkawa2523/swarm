"""Read and validate chemistry supplied for an exact set of collision data."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from hashlib import sha256
import math
from pathlib import Path
import re
from typing import Iterable

import yaml

from electron_swarm.core.yaml_input import _UniqueKeyLoader


CATALOG_SCHEMA = "swarm.electron_impact_catalog.v1"
ChannelKey = tuple[str, str, str, float | None]
Stoichiometry = tuple[tuple[str, int], ...]

# Numerical collision family, outgoing electrons, minimum heavy products.
# The non-dissociative categories require exactly one heavy product.
_CATEGORY_RULES = {
    "excitation": ("excitation", 1, 1),
    "superelastic": ("superelastic", 1, 1),
    "dissociation": ("excitation", 1, 2),
    "ionization": ("ionization", 2, 1),
    "dissociative_ionization": ("ionization", 2, 2),
    "attachment": ("attachment", 0, 1),
    "dissociative_attachment": ("attachment", 0, 2),
}
CHEMICAL_RATE_TYPES = frozenset(rule[0] for rule in _CATEGORY_RULES.values())


@dataclass(frozen=True, slots=True)
class Species:
    atoms: Stoichiometry
    charge: int


@dataclass(frozen=True, slots=True)
class Reaction:
    id: str
    key: ChannelKey
    category: str
    products: Stoichiometry
    resolution: str
    reference: str

    @property
    def reactants(self) -> Stoichiometry:
        return (("e", 1), (self.key[0], 1))

    def net(self, species: str) -> int:
        return dict(self.products).get(species, 0) - dict(self.reactants).get(
            species, 0
        )


@dataclass(frozen=True, slots=True)
class Dataset:
    name: str
    sha256: str
    reference: str


@dataclass(frozen=True, slots=True)
class ReactionCatalog:
    sha256: str
    datasets: tuple[Dataset, ...]
    species: dict[str, Species]
    reactions: dict[ChannelKey, Reaction]


@dataclass(frozen=True, slots=True)
class PathwayScreen:
    catalog_path: Path
    focus_species: str | None = None
    minimum_activity_s_inv: float = 0.0
    catalog_sha256: str = field(init=False)
    collision_file_sha256: tuple[str, ...] = field(init=False)

    def __post_init__(self) -> None:
        path = Path(self.catalog_path).resolve()
        catalog = read_reaction_catalog(path)
        if self.focus_species is not None and (
            not isinstance(self.focus_species, str)
            or self.focus_species not in catalog.species
        ):
            raise ValueError("reaction pathway focus_species must be a catalog species")
        minimum = _number(self.minimum_activity_s_inv, "minimum_activity_s_inv")
        object.__setattr__(self, "catalog_path", path)
        object.__setattr__(self, "minimum_activity_s_inv", minimum)
        object.__setattr__(self, "catalog_sha256", catalog.sha256)
        object.__setattr__(
            self,
            "collision_file_sha256",
            tuple(sorted(dataset.sha256 for dataset in catalog.datasets)),
        )

    def verify_collision_files(self, hashes: Iterable[str]) -> None:
        if tuple(sorted(hashes)) != self.collision_file_sha256:
            raise ValueError(
                "reaction catalog collision-file SHA-256 set differs from the Swarm inputs"
            )

    def plan(self) -> dict:
        return {
            "catalog": str(self.catalog_path),
            "catalog_sha256": self.catalog_sha256,
            "collision_file_sha256": list(self.collision_file_sha256),
            "focus_species": self.focus_species,
            "minimum_activity_s_inv": self.minimum_activity_s_inv,
            "scope": "direct_electron_impact_at_saved_dc_points",
            "channel_binding": "exact_match_to_native_rate_inventory",
            "grouping": "category"
            if self.focus_species is None
            else "focus_production_or_loss",
        }


def _mapping(value: object, label: str, fields: set[str] | None = None) -> dict:
    if not isinstance(value, dict) or (fields is not None and set(value) != fields):
        expected = f" with fields {sorted(fields)}" if fields is not None else ""
        raise ValueError(f"{label} must be a mapping{expected}")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError(
            f"{label} must be nonempty text without surrounding whitespace"
        )
    return value


def _number(value: object, label: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
        raise ValueError(f"{label} must be finite and nonnegative")
    return float(value)


def _counts(value: object, label: str) -> Stoichiometry:
    data = _mapping(value, label)
    if not data or any(type(count) is not int or count <= 0 for count in data.values()):
        raise ValueError(f"{label} must give positive integer counts")
    return tuple(sorted(data.items()))


def _read_species(raw: object) -> dict[str, Species]:
    species = {"e": Species((), -1)}
    for name, value in _mapping(raw, "species").items():
        _text(name, "species id")
        if name == "e":
            raise ValueError("species id e is reserved for the electron")
        value = _mapping(value, f"species {name}", {"atoms", "charge"})
        atoms = _counts(value["atoms"], f"species {name} atoms")
        if any(re.fullmatch(r"[A-Z][a-z]?", atom) is None for atom, _ in atoms):
            raise ValueError(f"species {name} atoms must use element symbols")
        charge = value["charge"]
        if type(charge) is not int:
            raise ValueError(f"species {name} charge must be an integer")
        species[name] = Species(atoms, charge)
    return species


def _balance(terms: Stoichiometry, species: dict[str, Species]) -> tuple[Counter, int]:
    atoms: Counter = Counter()
    charge = 0
    for name, count in terms:
        for element, quantity in species[name].atoms:
            atoms[element] += count * quantity
        charge += count * species[name].charge
    return atoms, charge


def _read_reaction(name: str, raw: object, species: dict[str, Species]) -> Reaction:
    label = f"reaction {name}"
    _text(name, "reaction id")
    data = _mapping(
        raw,
        label,
        {
            "target",
            "process",
            "threshold_eV",
            "category",
            "products",
            "resolution",
            "reference",
        },
    )
    target = _text(data["target"], f"{label} target")
    process = _text(data["process"], f"{label} process")
    category = _text(data["category"], f"{label} category")
    if category not in _CATEGORY_RULES:
        raise ValueError(
            f"{label} has unsupported electron-impact category {category!r}"
        )
    native, electrons, heavy = _CATEGORY_RULES[category]
    threshold = data["threshold_eV"]
    if threshold is not None or native != "attachment":
        threshold = _number(threshold, f"{label} threshold_eV")
    products = _counts(data["products"], f"{label} products")
    if (
        target not in species
        or target == "e"
        or any(item not in species for item, _ in products)
    ):
        raise ValueError(
            f"{label} requires declared heavy target and declared products"
        )
    count_e = dict(products).get("e", 0)
    count_heavy = sum(count for _, count in products) - count_e
    if count_e != electrons or count_heavy < heavy or (heavy == 1 and count_heavy != 1):
        raise ValueError(f"{label} product counts disagree with category {category}")
    resolution = _text(data["resolution"], f"{label} resolution")
    if resolution not in {"resolved", "lumped"}:
        raise ValueError(f"{label} resolution must be resolved or lumped")
    reaction = Reaction(
        name,
        (target, process, native, threshold),
        category,
        products,
        resolution,
        _text(data["reference"], f"{label} reference"),
    )
    if _balance(reaction.reactants, species) != _balance(products, species):
        raise ValueError(f"{label} violates element or charge balance")
    if dict(reaction.reactants) == dict(products):
        raise ValueError(f"{label} must change the target state or species")
    return reaction


def read_reaction_catalog(path: Path) -> ReactionCatalog:
    """Read one immutable YAML snapshot and verify its structural chemistry."""
    content = path.read_bytes()
    try:
        raw = yaml.load(content.decode("utf-8"), Loader=_UniqueKeyLoader)
    except (yaml.YAMLError, ValueError) as exc:
        raise ValueError(f"invalid reaction catalog {path}: {exc}") from exc
    data = _mapping(
        raw, "reaction catalog", {"schema", "datasets", "species", "reactions"}
    )
    if data["schema"] != CATALOG_SCHEMA:
        raise ValueError(f"reaction catalog schema must be {CATALOG_SCHEMA}")
    datasets = []
    for name, value in _mapping(data["datasets"], "datasets").items():
        _text(name, "dataset id")
        value = _mapping(value, f"dataset {name}", {"sha256", "reference"})
        digest = _text(value["sha256"], f"dataset {name} sha256")
        if re.fullmatch(r"[a-f0-9]{64}", digest) is None:
            raise ValueError(f"dataset {name} sha256 must be a lowercase SHA-256")
        datasets.append(
            Dataset(
                name, digest, _text(value["reference"], f"dataset {name} reference")
            )
        )
    if not datasets or len({item.sha256 for item in datasets}) != len(datasets):
        raise ValueError("reaction catalog requires distinct collision dataset hashes")
    species = _read_species(data["species"])
    reactions = {}
    for name, value in _mapping(data["reactions"], "reactions").items():
        reaction = _read_reaction(name, value, species)
        if reaction.key in reactions:
            raise ValueError(f"reaction {name} duplicates a collision channel")
        reactions[reaction.key] = reaction
    if not reactions:
        raise ValueError("reaction catalog must identify at least one channel")
    return ReactionCatalog(
        sha256(content).hexdigest(), tuple(datasets), species, reactions
    )
