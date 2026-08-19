import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from molgen.docking import parse_site_residues, binding_site_center, prepare_receptor_pdbqt, prepare_ligand_pdbqt
from molgen.docking_prep import embed_3d

REFERENCE_PDB = Path(__file__).resolve().parents[1] / "data" / "receptor" / "1M17.pdb"


def _skip_if_no_receptor():
    if not REFERENCE_PDB.exists():
        pytest.skip("data/receptor/1M17.pdb absent — lancer scripts/prepare_receptor.py d'abord")


def test_parse_site_residues_excludes_water():
    pdb_text = "SITE     1 AC1 14 HOH A  10  LEU A 694  ALA A 719  LEU A 764                    \n"
    residues = parse_site_residues(pdb_text)
    assert ("A", 694) in residues
    assert ("A", 719) in residues
    assert ("A", 764) in residues
    assert not any(r[1] == 10 for r in residues)  # HOH exclue


def test_binding_site_center_on_real_structure():
    _skip_if_no_receptor()
    pdb_text = REFERENCE_PDB.read_text()
    center = binding_site_center(pdb_text)
    assert center is not None
    x, y, z = center
    # sanity check large : les coordonnées doivent être dans une plage
    # physiquement raisonnable pour une structure cristallographique (pas
    # (0,0,0) ni des valeurs aberrantes issues d'un mauvais parsing)
    assert all(-200 < c < 200 for c in center)


def test_prepare_receptor_pdbqt_on_real_structure():
    _skip_if_no_receptor()
    out = Path("/tmp/test_receptor_pytest.pdbqt")
    ok = prepare_receptor_pdbqt(REFERENCE_PDB, out)
    assert ok
    assert out.exists()
    assert out.read_text().strip() != ""


def test_prepare_ligand_pdbqt():
    mol = embed_3d("C#Cc1ccc(Nc2ncnc3cc(OC)c(OCCCN4CCN(C)CC4)cc23)c(F)c1", mol_id="test")
    pdbqt = prepare_ligand_pdbqt(mol)
    assert pdbqt is not None
    assert "ATOM" in pdbqt
    assert "ROOT" in pdbqt  # bloc flexible Meeko/AutoDock


def test_dock_end_to_end():
    """Test d'intégration complet : parsing + préparation + docking réel,
    sur la vraie structure 1M17 si elle a été préparée localement."""
    _skip_if_no_receptor()
    from molgen.docking import dock

    receptor_pdbqt = Path("/tmp/test_receptor_dock.pdbqt")
    ok = prepare_receptor_pdbqt(REFERENCE_PDB, receptor_pdbqt)
    assert ok

    center = binding_site_center(REFERENCE_PDB.read_text())
    mol = embed_3d("COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCN1CCS(=O)(=O)CC1", mol_id="test")
    ligand_pdbqt = prepare_ligand_pdbqt(mol)

    score = dock(receptor_pdbqt, ligand_pdbqt, center, exhaustiveness=2)
    assert score is not None
    assert -20 < score < 0  # une affinité prédite est toujours négative, et rarement au-delà de -20 kcal/mol
