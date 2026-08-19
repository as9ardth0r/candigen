"""candigen — pipeline de génération et de criblage de molécules cibles (CADD)."""

from .properties import MoleculeRecord, compute_batch, compute_descriptors, load_smiles
from .filters import TPPProfile, enrich_and_filter

__all__ = [
    "MoleculeRecord",
    "compute_batch",
    "compute_descriptors",
    "load_smiles",
    "TPPProfile",
    "enrich_and_filter",
]

__version__ = "0.1.0"
