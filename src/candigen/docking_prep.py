"""
candigen.docking_prep
====================
Préparation des ligands pour le docking moléculaire : génération de
conformères 3D (ETKDGv3), optimisation MMFF94, export SDF/PDB, et
un gabarit de configuration AutoDock Vina.

Ce module s'arrête à l'export SDF/PDB : la conversion en .pdbqt (ajout des
charges Gasteiger, atomes rotatables) nécessite un outil externe non
inclus ici (Open Babel `obabel -ipdb -opdbqt` ou Meeko `mk_prepare_ligand.py`).
Voir le README pour l'installation de ces dépendances optionnelles.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def embed_3d(smiles: str, mol_id: str = "mol", seed: int = 42) -> Chem.Mol | None:
    """Génère un conformère 3D et l'optimise avec le champ de force MMFF94."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cid = AllChem.EmbedMolecule(mol, params)
    if cid < 0:
        # Repli si l'embedding échoue (molécules très contraintes)
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mol, params)
        if cid < 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    mol.SetProp("_Name", mol_id)
    return mol


def export_sdf(mol: Chem.Mol, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def mol_to_sdf_block(mol: Chem.Mol) -> str:
    """
    Retourne un bloc SDF valide (utile pour l'embarquer dans le JSON du site).

    IMPORTANT : `Chem.MolToMolBlock` retourne un bloc MOL (format molfile),
    pas un enregistrement SDF complet — il lui manque le terminateur `$$$$`
    requis par la spécification SDF. Sans lui, la plupart des parseurs SDF
    (dont celui de 3Dmol.js) ne trouvent aucun atome à charger, souvent sans
    lever d'erreur explicite (échec silencieux).
    """
    molblock = Chem.MolToMolBlock(mol)
    if not molblock.endswith("\n"):
        molblock += "\n"
    return molblock + "$$$$\n"


def export_pdb(mol: Chem.Mol, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToPDBFile(mol, str(path))


VINA_CONFIG_TEMPLATE = """\
receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}

center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

size_x = 20
size_y = 20
size_z = 20

exhaustiveness = 16
num_modes = 9
energy_range = 3
seed = 42
"""


def write_vina_config(
    path: str | Path,
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    center: tuple[float, float, float],
) -> None:
    """
    Écrit un fichier de config AutoDock Vina. Les coordonnées du site actif
    de l'EGFR (domaine kinase, poche ATP) doivent être déterminées à partir
    d'une structure cristallographique de référence (ex: PDB 1M17, 4HJO,
    ou une structure avec la mutation d'intérêt), pas codées en dur ici.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    content = VINA_CONFIG_TEMPLATE.format(
        receptor_pdbqt=receptor_pdbqt,
        ligand_pdbqt=ligand_pdbqt,
        center_x=center[0],
        center_y=center[1],
        center_z=center[2],
    )
    path.write_text(content, encoding="utf-8")
