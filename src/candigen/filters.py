"""
candigen.filters
==============
Filtres de criblage : règle de Lipinski (Ro5), profil cible (TPP) spécifique
à la cible, alertes structurales PAINS/Brenk, et score de synthétisabilité
(SA score, Ertl & Schuffenhauer 2009).

Le SA score utilise le script officiel RDKit Contrib (`vendor/sascorer.py`
+ `vendor/fpscores.pkl.gz`, cf. `scripts/fetch_vendor.sh`). S'il est absent,
on retombe sur une heuristique simplifiée basée sur la complexité
structurale (nombre de cycles, stéréocentres, hétéroatomes) — moins précise
mais qui ne bloque pas le pipeline.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from .properties import MoleculeRecord

VENDOR_DIR = Path(__file__).resolve().parents[2] / "vendor"

# --- SA score : script officiel si disponible, sinon heuristique de repli ---
_sascorer = None
if (VENDOR_DIR / "sascorer.py").exists() and (VENDOR_DIR / "fpscores.pkl.gz").exists():
    sys.path.insert(0, str(VENDOR_DIR))
    try:
        import sascorer as _sascorer  # type: ignore
    except Exception:
        _sascorer = None


def sa_score(mol: Chem.Mol) -> float:
    """Score de synthétisabilité, échelle 1 (facile) - 10 (très difficile)."""
    if _sascorer is not None:
        return round(_sascorer.calculateScore(mol), 2)
    # Heuristique de repli (approximative) : pénalise cycles fusionnés,
    # stéréocentres et hétéroatomes rares.
    ri = mol.GetRingInfo()
    n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    n_bridge = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
    score = 1.5 + 0.15 * ri.NumRings() + 0.3 * n_spiro + 0.3 * n_bridge + 0.2 * n_stereo
    return round(min(score, 10.0), 2)


# --- Catalogues d'alertes structurales ---
# PAINS : faux positifs de criblage (interfèrent avec le test, pas
# forcément dangereux en soi).
_pains_params = FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
_pains_catalog = FilterCatalog(_pains_params)

# BRENK : substructures réactives / toxicophores reconnus en chimie
# médicinale (groupes électrophiles, motifs génotoxiques ou instables...) —
# Brenk et al., ChemMedChem 2008. C'est le signal pertinent pour repérer
# une nocivité potentielle, à ne pas confondre avec PAINS.
_brenk_params = FilterCatalogParams()
_brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
_brenk_catalog = FilterCatalog(_brenk_params)


def structural_alerts(mol: Chem.Mol) -> list[str]:
    """Alertes PAINS (faux positifs de criblage)."""
    return [entry.GetDescription() for entry in _pains_catalog.GetMatches(mol)]


def toxicity_alerts(mol: Chem.Mol) -> list[str]:
    """Alertes BRENK (groupes réactifs / toxicophores connus) — signal de nocivité potentielle."""
    return [entry.GetDescription() for entry in _brenk_catalog.GetMatches(mol)]


def lipinski_violations(mol: Chem.Mol) -> int:
    """Règle des 5 de Lipinski : nombre de critères violés (0-1 = acceptable)."""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    return sum([mw > 500, logp > 5, hbd > 5, hba > 10])


@dataclass
class TPPProfile:
    """
    Target Product Profile — bornes physico-chimiques et ADMET utilisées
    pour cribler des candidats. Les valeurs par défaut correspondent à un
    inhibiteur ATP-compétitif de kinase à petite molécule, administré per os,
    calibré sur l'EGFR (poche à charnière type quinazoline/aminopyrimidine).
    """

    name: str = "EGFR kinase inhibitor (oral, Type I/covalent)"
    mw_range: tuple[float, float] = (250.0, 500.0)
    logp_range: tuple[float, float] = (1.0, 4.5)
    tpsa_range: tuple[float, float] = (40.0, 100.0)
    hbd_max: int = 3
    hba_max: int = 9
    rotb_max: int = 10
    sa_score_max: float = 4.0
    lipinski_violations_max: int = 1
    qed_min: float = 0.3
    forbid_pains: bool = True
    # BRENK reste INFORMATIF par défaut (pas de rejet automatique) : il
    # signale aussi des choix de conception légitimes et courants (ex.
    # warhead acrylamide d'un inhibiteur covalent, alcyne terminal...) —
    # cf. M2/M5 dans les candidats curés. La valeur est toujours calculée
    # et exposée (record.toxicity_alerts) pour une revue manuelle au cas
    # par cas. Passez à True pour un criblage automatique plus strict.
    forbid_toxicity_alerts: bool = False

    def evaluate(self, record: MoleculeRecord) -> tuple[bool, list[str]]:
        """Retourne (conforme: bool, raisons_de_rejet: list[str])."""
        reasons = []
        if not (self.mw_range[0] <= record.mw <= self.mw_range[1]):
            reasons.append(f"MW {record.mw} hors [{self.mw_range[0]}, {self.mw_range[1]}]")
        if not (self.logp_range[0] <= record.logp <= self.logp_range[1]):
            reasons.append(f"LogP {record.logp} hors [{self.logp_range[0]}, {self.logp_range[1]}]")
        if not (self.tpsa_range[0] <= record.tpsa <= self.tpsa_range[1]):
            reasons.append(f"TPSA {record.tpsa} hors [{self.tpsa_range[0]}, {self.tpsa_range[1]}]")
        if record.hbd > self.hbd_max:
            reasons.append(f"HBD {record.hbd} > {self.hbd_max}")
        if record.hba > self.hba_max:
            reasons.append(f"HBA {record.hba} > {self.hba_max}")
        if record.rotatable_bonds > self.rotb_max:
            reasons.append(f"RotB {record.rotatable_bonds} > {self.rotb_max}")
        if record.sa_score is not None and record.sa_score > self.sa_score_max:
            reasons.append(f"SA score {record.sa_score} > {self.sa_score_max}")
        if record.lipinski_violations is not None and record.lipinski_violations > self.lipinski_violations_max:
            reasons.append(f"Lipinski violations {record.lipinski_violations} > {self.lipinski_violations_max}")
        if record.qed < self.qed_min:
            reasons.append(f"QED {record.qed} < {self.qed_min}")
        if self.forbid_pains and record.pains_alerts:
            reasons.append(f"Alertes PAINS: {', '.join(record.pains_alerts)}")
        if self.forbid_toxicity_alerts and record.toxicity_alerts:
            reasons.append(f"Alertes BRENK (nocivité potentielle): {', '.join(record.toxicity_alerts)}")
        return (len(reasons) == 0, reasons)


def enrich_and_filter(records: list[MoleculeRecord], profile: Optional[TPPProfile] = None) -> list[MoleculeRecord]:
    """Calcule SA score / alertes / violations Lipinski et évalue chaque record contre le TPP."""
    profile = profile or TPPProfile()
    for r in records:
        mol = Chem.MolFromSmiles(r.canonical_smiles)
        r.sa_score = sa_score(mol)
        r.pains_alerts = structural_alerts(mol)
        r.toxicity_alerts = toxicity_alerts(mol)
        r.lipinski_violations = lipinski_violations(mol)
        passed, reasons = profile.evaluate(r)
        r.tpp_pass = passed
        r.notes = "; ".join(reasons) if reasons else "Conforme au TPP"
    return records
