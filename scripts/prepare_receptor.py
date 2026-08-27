#!/usr/bin/env python3
"""Prépare les structures réceptrices pour le docking (docking.py) :
télécharge le PDB depuis RCSB, repère le ligand co-cristallisé (pour
dériver automatiquement la boîte de recherche autour de sa position
réelle plutôt que de deviner des coordonnées), puis appelle
`mk_prepare_receptor.py` (meeko) pour produire le PDBQT et la
configuration de boîte Vina.

Nécessite un accès réseau vers files.rcsb.org, absent de la liste
blanche du sandbox utilisé pour construire ce dépôt (limitée à
PyPI/GitHub/npm) — ce script n'a donc pas pu être exécuté de bout en
bout ici. La logique d'extraction du ligand (receptor_prep.py) est
testée séparément sans réseau — voir tests/test_receptor_prep.py.
Fonctionnera normalement en local ou dans explore.yml (les runners
GitHub Actions ont un accès réseau complet, contrairement à ce sandbox).

Les identifiants PDB ci-dessous sont des points de départ à vérifier
(résolution, complétude du domaine kinase, présence d'un ligand
co-cristallisé pertinent) avant tout usage réel — pas un choix validé.

Usage prévu :
    python scripts/prepare_receptor.py --target EGFR_WT
    python scripts/prepare_receptor.py --target all
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from candigen.receptor_prep import (
    extract_ligand_pdb,
    find_cocrystallized_ligand,
    find_all_residue_instances,
    format_delete_residues,
)  # noqa: E402

CANDIDATE_STRUCTURES = {
    "EGFR_WT": "1M17",              # domaine kinase EGFR + erlotinib
    "EGFR_T790M": "2JIV",           # mutant de résistance T790M
    "EGFR_T790M_C797S": "5D41",     # double mutant, cible de l'osimertinib
    "HER2": "3PP0",                 # domaine kinase HER2/ErbB2
}

REPO_ROOT = Path(__file__).resolve().parent.parent
RECEPTOR_DIR = REPO_ROOT / "data" / "receptors"


def fetch_pdb(pdb_id: str, out_path: Path) -> None:
    import requests
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    out_path.write_text(response.text)


def prepare_receptor(pdb_id: str, name: str, padding: float, allow_bad_res: bool, default_altloc: str) -> None:
    raw_path = RECEPTOR_DIR / f"{name.lower()}_raw.pdb"
    print(f"[prepare_receptor] téléchargement {pdb_id} -> {raw_path}")
    fetch_pdb(pdb_id, raw_path)

    ligand_hit = find_cocrystallized_ligand(raw_path)
    output_basename = str(RECEPTOR_DIR / name.lower())
    cmd = [
        "mk_prepare_receptor.py",
        "--read_pdb", str(raw_path),
        "-o", output_basename,
        "-p",  # écrit le PDBQT
        "-v",  # écrit la config de boîte Vina
        "--padding", str(padding),
        # Les structures cristallographiques réelles ont souvent des résidus
        # à localisation alternative (ex. 1M17 : A:751, A:831) — deux
        # conformations modélisées pour la même chaîne latérale, ambiguïté
        # que meeko refuse de trancher seul. 'A' est la convention PDB pour
        # la conformation principale (occupation la plus élevée), donc un
        # défaut raisonnable plutôt que de supprimer ces résidus.
        "--default_altloc", default_altloc,
    ]
    if allow_bad_res:
        cmd.append("-a")

    if ligand_hit is not None:
        chain_name, res_name, res_seq, n_atoms = ligand_hit
        print(f"[prepare_receptor] ligand co-cristallisé repéré : "
              f"{res_name} (chaîne {chain_name}, résidu {res_seq}, {n_atoms} atomes)")
        ligand_path = RECEPTOR_DIR / f"{name.lower()}_ref_ligand.pdb"
        extract_ligand_pdb(raw_path, chain_name, res_name, res_seq, ligand_path)
        cmd += ["--box_enveloping", str(ligand_path)]

        # Le récepteur préparé doit être exempt du ligand co-cristallisé —
        # sinon il reste "en dur" dans la structure rigide et bloque
        # physiquement sa propre poche de liaison pour tout futur docking.
        # On l'exclut systématiquement (pas seulement quand meeko échoue à
        # lui générer un template, ex. HKI dans 2JIV) et sur TOUTES ses
        # copies (une structure peut avoir plusieurs chaînes protéiques
        # dans l'unité asymétrique, chacune avec sa propre copie du ligand).
        instances = find_all_residue_instances(raw_path, res_name)
        delete_arg = format_delete_residues(instances)
        cmd += ["--delete_residues", delete_arg]
        print(f"[prepare_receptor] exclusion du récepteur : {res_name} "
              f"({len(instances)} copie(s) : {delete_arg})")
    else:
        print(f"[prepare_receptor] ATTENTION — aucun ligand co-cristallisé trouvé pour {name} : "
              f"la boîte de recherche doit être définie manuellement (--box_center / --box_size), "
              f"pas dérivée automatiquement.")

    print(f"[prepare_receptor] exécution : {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"[prepare_receptor] {name} prêt : "
          f"{output_basename}.pdbqt + {output_basename}_box.txt")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", choices=list(CANDIDATE_STRUCTURES) + ["all"], default="all")
    parser.add_argument("--padding", type=float, default=4.0,
                         help="marge (Å) ajoutée autour du ligand de référence pour la boîte de recherche")
    parser.add_argument("--allow-bad-residues", action=argparse.BooleanOptionalAction, default=True,
                         help="supprime les résidus qui ne correspondent à aucun template plutôt que de "
                              "stopper sur erreur (voir -a de mk_prepare_receptor.py). Activé par défaut : "
                              "c'est la pratique recommandée par meeko lui-même pour des structures PDB "
                              "réelles, presque toujours partiellement résolues (--no-allow-bad-residues "
                              "pour désactiver).")
    parser.add_argument("--default-altloc", type=str, default="A",
                         help="conformation alternative par défaut pour les résidus qui en ont "
                              "plusieurs dans le cristal (voir --default_altloc de mk_prepare_receptor.py)")
    args = parser.parse_args()

    targets = list(CANDIDATE_STRUCTURES) if args.target == "all" else [args.target]
    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)

    # Chaque cible est préparée à partir d'une vraie structure cristallo-
    # graphique déposée par des équipes différentes, à des époques
    # différentes, avec ses propres défauts (boucles manquantes, résidus à
    # localisation alternative, chaînes multiples...). Il n'existe pas de
    # combinaison de flags qui absorbe tous les cas par avance — et le
    # reste du pipeline (dock_top_candidates.py) sait déjà gérer une cible
    # non préparée. Un échec sur UNE structure ne doit donc plus bloquer
    # les trois autres : on isole les échecs plutôt que de les empêcher.
    failures: list[str] = []
    for name in targets:
        pdb_id = CANDIDATE_STRUCTURES[name]
        try:
            prepare_receptor(pdb_id, name, args.padding, args.allow_bad_residues, args.default_altloc)
        except Exception as exc:
            print(f"[prepare_receptor] ÉCHEC pour {name} ({pdb_id}) : {exc}")
            print(f"[prepare_receptor] {name} ignoré — poursuite avec les cibles restantes.")
            failures.append(name)

    if failures:
        print(f"\n[prepare_receptor] {len(failures)}/{len(targets)} cible(s) en échec : {', '.join(failures)}")
        print("[prepare_receptor] Ces structures ont probablement besoin d'une inspection manuelle "
              "(boucles manquantes créant des liaisons parasites, résidu non standard sans template...) "
              "— pas d'un nouveau flag générique. dock_top_candidates.py ignore proprement une cible "
              "sans récepteur préparé.")
    if len(failures) == len(targets):
        sys.exit(1)  # tout a échoué : ça, c'est un vrai problème (réseau, RCSB down, etc.)


if __name__ == "__main__":
    main()
