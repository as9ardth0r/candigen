import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest
from candigen.properties import compute_descriptors, InvalidSMILESError
from candigen.filters import TPPProfile, lipinski_violations, sa_score
from rdkit import Chem


def test_aspirin_descriptors():
    r = compute_descriptors("aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O")
    assert r.formula == "C9H8O4"
    assert 179 < r.mw < 181
    assert r.hbd == 1
    assert r.hba == 3  # définition RDKit (NOCount) — diffère d'un simple compte O+N


def test_invalid_smiles_raises():
    with pytest.raises(InvalidSMILESError):
        compute_descriptors("bad", "not_a_smiles(((")


def test_lipinski_violations_large_molecule():
    # Une molécule volumineuse et très lipophile doit violer plusieurs règles
    huge_lipophilic = "C" * 40
    mol = Chem.MolFromSmiles(huge_lipophilic)
    assert lipinski_violations(mol) >= 1


def test_sa_score_in_range():
    mol = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    score = sa_score(mol)
    assert 1.0 <= score <= 10.0


def test_tpp_profile_evaluates():
    r = compute_descriptors("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    profile = TPPProfile()
    passed, reasons = profile.evaluate(r)
    assert isinstance(passed, bool)
    assert isinstance(reasons, list)


def test_generator_produces_valid_analogs():
    from candigen.generator import enumerate_analogs, ANILINE_HINGE_SUBSTITUENTS, SOLUBILIZING_ARMS

    scaffold = "c1ccc2ncnc(N[*:1])c2c1[*:2]"
    analogs = enumerate_analogs(scaffold, {1: ANILINE_HINGE_SUBSTITUENTS, 2: SOLUBILIZING_ARMS})
    assert len(analogs) == len(ANILINE_HINGE_SUBSTITUENTS.fragments) * len(SOLUBILIZING_ARMS.fragments)
    for mol_id, smi in analogs:
        assert Chem.MolFromSmiles(smi) is not None, f"{mol_id} produced invalid SMILES"


def test_full_library_is_large_and_valid():
    from candigen.generator import enumerate_library, SCAFFOLD_LIBRARY, ANILINE_HINGE_SUBSTITUENTS, SOLUBILIZING_ARMS

    library = enumerate_library(SCAFFOLD_LIBRARY, {1: ANILINE_HINGE_SUBSTITUENTS, 2: SOLUBILIZING_ARMS})
    expected = len(SCAFFOLD_LIBRARY) * len(ANILINE_HINGE_SUBSTITUENTS.fragments) * len(SOLUBILIZING_ARMS.fragments)
    assert len(library) == expected
    assert expected >= 100  # une "vraie" bibliothèque, pas une poignée de molécules
    invalid = [mol_id for mol_id, smi in library if Chem.MolFromSmiles(smi) is None]
    assert not invalid, f"SMILES invalides générés : {invalid}"
