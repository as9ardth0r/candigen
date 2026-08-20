"""
candigen.properties
=================
Chargement de SMILES et calcul des descripteurs physico-chimiques
(MW, LogP, TPSA, HBD, HBA, rotatable bonds, QED) via RDKit.

C'est le module "cœur de calcul" du pipeline : toute la suite
(filtres, export, site web) consomme des `MoleculeRecord`.
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Iterable, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors

# RDKit est bruyant sur les SMILES invalides -> on gère nous-mêmes le logging
RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


@dataclass
class MoleculeRecord:
    """Représentation d'une molécule et de ses descripteurs calculés."""

    id: str
    smiles: str
    canonical_smiles: str
    mw: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    num_rings: int
    num_aromatic_rings: int
    qed: float
    formula: str
    sa_score: Optional[float] = None
    lipinski_violations: Optional[int] = None
    pains_alerts: list[str] = field(default_factory=list)
    toxicity_alerts: list[str] = field(default_factory=list)  # alertes BRENK (réactivité/nocivité potentielle)
    tpp_pass: Optional[bool] = None
    notes: str = ""
    source: str = "generated"
    recipe: Optional[dict] = None  # {"scaffold":..., "aniline":..., "solubilizer":...} pour source="generated"
    first_seen: Optional[str] = None  # date ISO de première apparition dans le hall of fame
    fitness: Optional[float] = None  # score composite (voir candigen.evolve.fitness)
    docking_score: Optional[float] = None  # kcal/mol (AutoDock Vina), None si pas encore docké
    is_novel: Optional[bool] = None  # None = pas vérifié ; True = confirmé absent de PubChem/ChEMBL ; False = déjà connue
    pubchem_cid: Optional[int] = None
    chembl_id: Optional[str] = None
    chemical_name: Optional[str] = None  # nom IUPAC calculé par PubChem, si trouvé (cf. candigen.novelty) — None pour une molécule absente de PubChem (typiquement les molécules générées)
    retrosynthesis_route_found: Optional[bool] = None  # None = pas évaluée ; True/False = route trouvée ou non
    retrosynthesis_n_routes: Optional[int] = None  # nombre de routes retournées par AiZynthFinder

    def to_dict(self) -> dict:
        return asdict(self)


class InvalidSMILESError(ValueError):
    """Levée quand une chaîne SMILES ne peut pas être parsée par RDKit."""


def load_smiles(source: str | Path) -> list[tuple[str, str]]:
    """
    Charge une liste de (id, smiles) depuis :
      - un fichier .smi / .txt  (une entrée par ligne : "SMILES id")
      - un fichier .csv         (colonnes 'id' et 'smiles')
      - une simple chaîne SMILES unique

    Retourne une liste de tuples (id, smiles).
    """
    path = Path(source) if isinstance(source, (str, Path)) else None

    if path is not None and path.exists():
        if path.suffix == ".csv":
            with path.open(newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                return [(row["id"], row["smiles"]) for row in reader]
        else:  # .smi / .txt
            entries = []
            with path.open(encoding="utf-8") as f:
                for i, line in enumerate(f):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    smi = parts[0]
                    mol_id = parts[1] if len(parts) > 1 else f"mol_{i:03d}"
                    entries.append((mol_id, smi))
            return entries

    # Fallback : chaîne SMILES unique passée directement
    return [("mol_000", str(source))]


def compute_descriptors(mol_id: str, smiles: str) -> MoleculeRecord:
    """Calcule tous les descripteurs pour un SMILES donné. Lève InvalidSMILESError si invalide."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSMILESError(f"SMILES invalide pour '{mol_id}': {smiles!r}")

    return MoleculeRecord(
        id=mol_id,
        smiles=smiles,
        canonical_smiles=Chem.MolToSmiles(mol),
        mw=round(Descriptors.MolWt(mol), 2),
        logp=round(Crippen.MolLogP(mol), 2),
        tpsa=round(rdMolDescriptors.CalcTPSA(mol), 2),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
        num_rings=rdMolDescriptors.CalcNumRings(mol),
        num_aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
        qed=round(QED.qed(mol), 3),
        formula=rdMolDescriptors.CalcMolFormula(mol),
    )


def compute_batch(entries: Iterable[tuple[str, str]], skip_invalid: bool = True) -> list[MoleculeRecord]:
    """Calcule les descripteurs pour un lot de (id, smiles), en journalisant les échecs."""
    records = []
    for mol_id, smi in entries:
        try:
            records.append(compute_descriptors(mol_id, smi))
        except InvalidSMILESError as exc:
            logger.warning(str(exc))
            if not skip_invalid:
                raise
    return records


def export_json(records: list[MoleculeRecord], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump([r.to_dict() for r in records], f, indent=2, ensure_ascii=False)


def export_csv(records: list[MoleculeRecord], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        return
    fieldnames = list(records[0].to_dict().keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in records:
            row = r.to_dict()
            row["pains_alerts"] = ";".join(row["pains_alerts"])
            row["toxicity_alerts"] = ";".join(row["toxicity_alerts"])
            writer.writerow(row)
