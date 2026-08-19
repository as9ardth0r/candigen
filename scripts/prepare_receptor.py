#!/usr/bin/env python3
"""
Prépare la structure receptrice pour le docking — À LANCER UNE SEULE FOIS
(le résultat est committé dans data/receptor/ et réutilisé à chaque run du
pipeline, pas retéléchargé/reconverti systématiquement).

Structure utilisée : PDB 1M17 (EGFR, domaine kinase, avec erlotinib —
Stamos et al., J. Biol. Chem. 2002). Choisie parce que c'est la structure
de référence historique pour les inhibiteurs réversibles de Type I à cœur
quinazoline/aminopyrimidine — le chimiotype majoritaire de ce projet.
Pour cribler contre un mutant (T790M/L858R, résistance), remplacer
PDB_ID par "4HJO" ou une structure équivalente et relancer ce script.

Usage :
    python scripts/prepare_receptor.py
"""

from __future__ import annotations

import subprocess
import sys
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from candigen.docking import binding_site_center, prepare_receptor_pdbqt

PDB_ID = "1M17"
RECEPTOR_DIR = ROOT / "data" / "receptor"


def main() -> None:
    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)
    pdb_path = RECEPTOR_DIR / f"{PDB_ID}.pdb"

    print(f"[1/4] Téléchargement de {PDB_ID}.pdb depuis RCSB...")
    url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
    urllib.request.urlretrieve(url, pdb_path)
    pdb_text = pdb_path.read_text()
    print(f"      {len(pdb_text.splitlines())} lignes téléchargées")

    print("[2/4] Calcul du centre de la poche de liaison (records SITE)...")
    center = binding_site_center(pdb_text)
    if center is None:
        print("      ERREUR : aucun centre trouvé — vérifier les records SITE du PDB")
        sys.exit(1)
    print(f"      Centre : {center}")

    print("[3/4] Conversion du récepteur en PDBQT (Open Babel)...")
    receptor_pdbqt = RECEPTOR_DIR / f"{PDB_ID}.pdbqt"
    ok = prepare_receptor_pdbqt(pdb_path, receptor_pdbqt)
    if not ok:
        print("      ERREUR : conversion échouée — le paquet 'openbabel' (apt) est-il installé ?")
        sys.exit(1)
    print(f"      -> {receptor_pdbqt}")

    print("[4/4] Sauvegarde de la configuration...")
    config = {
        "pdb_id": PDB_ID,
        "center": list(center),
        "box_size": [24.0, 20.0, 24.0],
        "receptor_pdbqt": str(receptor_pdbqt.relative_to(ROOT)),
    }
    config_path = RECEPTOR_DIR / "config.json"
    config_path.write_text(json.dumps(config, indent=2))
    print(f"      -> {config_path}")
    print("\nTerminé. Committez data/receptor/ dans le dépôt pour que le pipeline l'utilise.")


if __name__ == "__main__":
    main()
