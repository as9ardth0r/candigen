"""
molgen.docking
==============
Docking moléculaire réel contre une structure receptrice EGFR (PDB), pour
compléter le criblage 2D (TPP/SA/PAINS/BRENK) par un score d'affinité basé
sur la géométrie 3D. Testé de bout en bout (voir README) avec la structure
1M17 (EGFR + erlotinib) :
  - extraction du centre de la poche depuis les records SITE du PDB,
  - préparation du récepteur en PDBQT via Open Babel,
  - préparation du ligand en PDBQT via Meeko,
  - docking avec AutoDock Vina (bindings Python officiels).

Le récepteur (fichier .pdbqt + centre de la poche) est préparé UNE FOIS par
`scripts/prepare_receptor.py` et committé dans data/receptor/ — le docking
proprement dit réutilise ces fichiers à chaque run, sans retélécharger ni
reconvertir le PDB à chaque fois.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from rdkit import Chem


def parse_site_residues(pdb_text: str) -> list[tuple[str, int]]:
    """
    Extrait les résidus listés dans les records SITE d'un fichier PDB
    (résidus en contact avec le ligand co-cristallisé), en excluant l'eau.
    Retourne une liste de (chain, resnum).
    """
    residues: list[tuple[str, int]] = []
    for line in pdb_text.splitlines():
        if not line.startswith("SITE"):
            continue
        for i in range(4):  # jusqu'à 4 résidus par ligne SITE
            start = 18 + i * 11
            chunk = line[start:start + 11]
            if len(chunk) < 10 or not chunk.strip():
                continue
            resname = chunk[0:3].strip()
            chain = chunk[4:5].strip()
            resnum_str = chunk[5:10].strip()
            if resname == "HOH" or not resnum_str:
                continue
            try:
                resnum = int(resnum_str)
            except ValueError:
                continue
            residues.append((chain, resnum))
    return residues


def binding_site_center(pdb_text: str) -> tuple[float, float, float] | None:
    """
    Calcule le centre géométrique (centroïde) de la poche de liaison à
    partir des résidus SITE et de leurs coordonnées ATOM dans le même
    fichier. Retourne None si aucun résidu SITE ou aucune coordonnée
    correspondante n'est trouvée.
    """
    site_residues = set(parse_site_residues(pdb_text))
    if not site_residues:
        return None

    xs, ys, zs = [], [], []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        chain = line[21:22].strip()
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        if (chain, resnum) not in site_residues:
            continue
        try:
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        except ValueError:
            continue

    if not xs:
        return None
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def prepare_receptor_pdbqt(pdb_path: str | Path, output_path: str | Path) -> bool:
    """
    Convertit un PDB de récepteur en PDBQT rigide via Open Babel
    (`obabel input.pdb -O output.pdbqt -xr`). Retourne True si la
    conversion a réussi. Nécessite le paquet `openbabel` (apt), pas
    seulement le module Python — voir requirements du workflow CI.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        result = subprocess.run(
            ["obabel", str(pdb_path), "-O", str(output_path), "-xr"],
            capture_output=True, text=True, timeout=120,
        )
        return result.returncode == 0 and output_path.exists()
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def prepare_ligand_pdbqt(mol: Chem.Mol) -> str | None:
    """Convertit un ligand (molécule RDKit avec conformère 3D) en chaîne PDBQT via Meeko."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    try:
        prep = MoleculePreparation()
        mol_setups = prep.prepare(mol)
        pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
        return pdbqt_string if is_ok else None
    except Exception:
        return None


def dock(
    receptor_pdbqt_path: str | Path,
    ligand_pdbqt_string: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float] = (24.0, 20.0, 24.0),
    exhaustiveness: int = 4,
) -> float | None:
    """
    Lance AutoDock Vina et retourne le meilleur score d'affinité prédit
    (kcal/mol — plus négatif = liaison plus favorable). Retourne None en
    cas d'échec (récepteur/ligand invalide, etc.) plutôt que de lever une
    exception, pour ne pas interrompre le criblage d'un lot entier à cause
    d'une molécule.

    Coût mesuré sur ce projet (runner 2 cœurs) : ~15s à exhaustiveness=2,
    ~30s à exhaustiveness=4, ~55s à exhaustiveness=8, à peu près linéaire.
    exhaustiveness=4 est un compromis raisonnable vitesse/qualité pour un
    run quotidien automatisé ; monter à 8-16 pour une analyse ponctuelle
    plus poussée sur une poignée de molécules.
    """
    from vina import Vina

    try:
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(str(receptor_pdbqt_path))
        v.set_ligand_from_string(ligand_pdbqt_string)
        v.compute_vina_maps(center=list(center), box_size=list(box_size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        scores = v.energies()
        if scores is None or len(scores) == 0:
            return None
        return round(float(scores[0][0]), 2)
    except Exception:
        return None
