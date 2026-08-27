import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest
from candigen.receptor_prep import (
    find_cocrystallized_ligand,
    extract_ligand_pdb,
    find_all_residue_instances,
    format_delete_residues,
)


# PDB minimal synthétique : eau + un ion (ZN) + un additif de cristallisation
# (GOL, glycérol, 6 atomes) + un "vrai" ligand (LG1, 12 atomes) sur la chaîne A.
# Pas de fichier PDB réel téléchargé — files.rcsb.org n'est pas dans la liste
# blanche réseau de ce sandbox (cf. docstring de scripts/prepare_receptor.py).
_SYNTHETIC_PDB = """\
HETATM  501  O   HOH A 301      10.000  10.000  10.000  1.00 20.00           O
HETATM  502  O   HOH A 302      11.000  10.000  10.000  1.00 20.00           O
HETATM  503 ZN    ZN A 401      15.000  15.000  15.000  1.00 15.00          ZN
HETATM  504  C1  GOL A 402      20.000  20.000  20.000  1.00 25.00           C
HETATM  505  C2  GOL A 402      20.500  20.500  20.500  1.00 25.00           C
HETATM  506  C3  GOL A 402      21.000  21.000  21.000  1.00 25.00           C
HETATM  507  O1  GOL A 402      21.500  21.500  21.500  1.00 25.00           O
HETATM  508  O2  GOL A 402      22.000  22.000  22.000  1.00 25.00           O
HETATM  509  O3  GOL A 402      22.500  22.500  22.500  1.00 25.00           O
HETATM  510  N1  LG1 A 501      30.000  30.000  30.000  1.00 18.00           N
HETATM  511  C1  LG1 A 501      30.500  30.500  30.500  1.00 18.00           C
HETATM  512  C2  LG1 A 501      31.000  31.000  31.000  1.00 18.00           C
HETATM  513  C3  LG1 A 501      31.500  31.500  31.500  1.00 18.00           C
HETATM  514  C4  LG1 A 501      32.000  32.000  32.000  1.00 18.00           C
HETATM  515  C5  LG1 A 501      32.500  32.500  32.500  1.00 18.00           C
HETATM  516  N2  LG1 A 501      33.000  33.000  33.000  1.00 18.00           N
HETATM  517  O1  LG1 A 501      33.500  33.500  33.500  1.00 18.00           O
HETATM  518  C6  LG1 A 501      34.000  34.000  34.000  1.00 18.00           C
HETATM  519  C7  LG1 A 501      34.500  34.500  34.500  1.00 18.00           C
HETATM  520  C8  LG1 A 501      35.000  35.000  35.000  1.00 18.00           C
HETATM  521  C9  LG1 A 501      35.500  35.500  35.500  1.00 18.00           C
"""


@pytest.fixture
def synthetic_pdb(tmp_path):
    p = tmp_path / "synthetic.pdb"
    p.write_text(_SYNTHETIC_PDB)
    return p


def test_finds_ligand_not_water_ion_or_cryoprotectant(synthetic_pdb):
    result = find_cocrystallized_ligand(synthetic_pdb)
    assert result == ("A", "LG1", "501", 12)


def test_returns_none_when_no_ligand_candidate(tmp_path):
    p = tmp_path / "only_water.pdb"
    p.write_text(
        "HETATM  501  O   HOH A 301      10.000  10.000  10.000  1.00 20.00           O\n"
        "HETATM  502 ZN    ZN A 401      15.000  15.000  15.000  1.00 15.00          ZN\n"
    )
    assert find_cocrystallized_ligand(p) is None


def test_extract_ligand_pdb_writes_only_matching_atoms(synthetic_pdb, tmp_path):
    out = tmp_path / "ligand.pdb"
    extract_ligand_pdb(synthetic_pdb, "A", "LG1", "501", out)
    content = out.read_text()
    lines = [ln for ln in content.splitlines() if ln.startswith("HETATM")]
    assert len(lines) == 12
    assert all("LG1" in ln for ln in lines)
    assert content.rstrip().endswith("END")


def test_extract_ligand_pdb_raises_on_no_match(synthetic_pdb, tmp_path):
    out = tmp_path / "ligand.pdb"
    with pytest.raises(ValueError):
        extract_ligand_pdb(synthetic_pdb, "A", "NOPE", "999", out)


# PDB synthétique avec le MÊME ligand sur deux chaînes différentes — simule
# 2JIV, où HKI apparaît à la fois sur la chaîne A et la chaîne B (deux
# copies de la protéine dans l'unité asymétrique).
_TWO_CHAIN_PDB = """\
HETATM  601  N1  LIG A 501      10.000  10.000  10.000  1.00 20.00           N
HETATM  602  C1  LIG A 501      10.500  10.500  10.500  1.00 20.00           C
HETATM  603  N1  LIG B 501      50.000  50.000  50.000  1.00 20.00           N
HETATM  604  C1  LIG B 501      50.500  50.500  50.500  1.00 20.00           C
"""


def test_find_all_residue_instances_across_chains(tmp_path):
    p = tmp_path / "two_chain.pdb"
    p.write_text(_TWO_CHAIN_PDB)
    instances = find_all_residue_instances(p, "LIG")
    assert set(instances) == {("A", "501"), ("B", "501")}


def test_format_delete_residues_groups_by_chain():
    result = format_delete_residues([("A", "999")])
    assert result == "A:999"

    result2 = format_delete_residues([("A", "501"), ("B", "501")])
    assert result2 == "A:501,B:501"

    result3 = format_delete_residues([("A", "15"), ("A", "16"), ("B", "42")])
    assert result3 == "A:15,16,B:42"
