"""
candigen.receptor_prep
========================
Repérage et extraction du ligand co-cristallisé dans un fichier PDB brut
(téléchargé depuis RCSB) : sert à dériver automatiquement la boîte de
recherche AutoDock Vina autour de sa position réelle plutôt que de deviner
des coordonnées à la main.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

# Résidus HETATM à ignorer car ce ne sont jamais le ligand d'intérêt pour
# l'arrimage moléculaire : eau, ions courants, additifs de cristallisation
# / cryoprotection, et résidus standards (au cas où un résidu modifié
# serait marqué HETATM). Liste pragmatique, pas exhaustive — à étendre si
# un nouveau faux positif apparaît sur une structure donnée.
_IGNORED_RESNAMES = {
    # eau
    "HOH", "WAT", "DOD",
    # ions courants
    "NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE", "FE2", "CO", "NI",
    "CU", "CU1", "CD", "HG", "BA", "CS", "LI", "RB", "SR", "AL", "PB",
    "AG", "PT", "AU", "GD", "SM", "YB", "LA", "CE", "TB", "IOD", "BR",
    # additifs de cristallisation / cryoprotection courants
    "GOL", "EDO", "PEG", "PG4", "1PE", "2PE", "MPD", "DMS", "SO4",
    "PO4", "ACT", "TRS", "BME", "MRD", "IPA", "FMT", "CIT", "IMD",
    "EPE", "BTB", "BCT", "NO3", "ACY", "CAC", "MES", "BEZ", "PGE",
    "P6G", "1PG", "DIO", "TAM", "UNK",
    # acides aminés et nucléotides standards
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT",
}


def _read_hetatm_lines(pdb_path: str | Path) -> list[str]:
    lines = Path(pdb_path).read_text(errors="replace").splitlines()
    return [ln for ln in lines if ln.startswith("HETATM")]


def find_cocrystallized_ligand(pdb_path: str | Path) -> tuple[str, str, str, int] | None:
    """Scanne un fichier PDB brut et retourne le ligand co-cristallisé le
    plus probable, sous la forme (chain_id, res_name, res_seq, n_atoms).

    Heuristique : parmi tous les groupes HETATM distincts (regroupés par
    chaîne + nom de résidu + numéro de résidu) qui ne sont ni de l'eau, ni
    un ion courant, ni un additif de cristallisation connu (voir
    `_IGNORED_RESNAMES`), on retient celui qui a le plus d'atomes — un vrai
    ligand pharmacologique (ex. un inhibiteur de kinase) a typiquement
    15 à 40 atomes lourds, largement plus que les additifs filtrés.

    Retourne None si aucun candidat plausible n'est trouvé (structure sans
    ligand co-cristallisé exploitable — la boîte de docking devra être
    définie manuellement via --box_center / --box_size).
    """
    groups: dict[tuple[str, str, str], int] = defaultdict(int)
    for line in _read_hetatm_lines(pdb_path):
        res_name = line[17:20].strip()
        chain_id = line[21].strip()
        res_seq = line[22:26].strip()
        if res_name in _IGNORED_RESNAMES:
            continue
        groups[(chain_id, res_name, res_seq)] += 1

    if not groups:
        return None

    (chain_id, res_name, res_seq), n_atoms = max(groups.items(), key=lambda kv: kv[1])
    return chain_id, res_name, res_seq, n_atoms


def extract_ligand_pdb(
    pdb_path: str | Path,
    chain_id: str,
    res_name: str,
    res_seq: str,
    output_path: str | Path,
) -> None:
    """Extrait les lignes HETATM d'un ligand précis (identifié via
    `find_cocrystallized_ligand`) dans un nouveau fichier PDB minimal,
    utilisable comme référence pour définir la boîte de recherche Vina
    (`--box_enveloping`)."""
    matching = [
        line for line in _read_hetatm_lines(pdb_path)
        if line[21].strip() == chain_id
        and line[17:20].strip() == res_name
        and line[22:26].strip() == res_seq
    ]

    if not matching:
        raise ValueError(
            f"Aucun atome trouvé pour {res_name} (chaîne {chain_id}, résidu {res_seq}) dans {pdb_path}"
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(matching) + "\nEND\n")


def find_all_residue_instances(pdb_path: str | Path, res_name: str) -> list[tuple[str, str]]:
    """Retourne (chain_id, res_seq) pour CHAQUE occurrence de `res_name`
    dans le fichier — un ligand co-cristallisé peut apparaître plusieurs
    fois (ex. plusieurs copies de la protéine dans l'unité asymétrique :
    2JIV a le ligand HKI à la fois sur la chaîne A et la chaîne B). Sert à
    s'assurer qu'AUCUNE copie ne reste dans le récepteur préparé — voir
    format_delete_residues.
    """
    instances: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for line in _read_hetatm_lines(pdb_path):
        if line[17:20].strip() != res_name:
            continue
        key = (line[21].strip(), line[22:26].strip())
        if key not in seen:
            seen.add(key)
            instances.append(key)
    return instances


def format_delete_residues(instances: list[tuple[str, str]]) -> str:
    """Formate une liste de (chain_id, res_seq) au format attendu par
    `mk_prepare_receptor.py --delete_residues` (ex. 'A:350,B:15,16,17' —
    un groupe chain:res[,res...] par chaîne, groupes séparés par des
    virgules, voir `mk_prepare_receptor.py --help`)."""
    by_chain: dict[str, list[str]] = {}
    for chain_id, res_seq in instances:
        by_chain.setdefault(chain_id, []).append(res_seq)
    return ",".join(f"{chain}:{','.join(resnums)}" for chain, resnums in by_chain.items())
