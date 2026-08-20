#!/usr/bin/env python3
"""
Prépare une ou plusieurs structures receptrices pour le docking multi-cible.
Idempotent par cible : relancer ce script ne re-télécharge/reconvertit QUE
les cibles absentes de data/receptor/<PDB_ID>/config.json — les cibles déjà
préparées et committées sont sautées automatiquement (pas besoin de logique
conditionnelle côté CI, contrairement à avant).

Pour ajouter/retirer une cible, éditer TARGETS ci-dessous et relancer.
Chaque cible est indépendante : data/receptor/<PDB_ID>/{<PDB_ID>.pdb,
<PDB_ID>.pdbqt, config.json}.

Usage :
    python scripts/prepare_receptor.py
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from candigen.docking import binding_site_center, prepare_receptor_pdbqt

# Une entrée par cible. target_name est la seule source de vérité pour le
# nom affiché dans le dashboard (candigen.export.read_target_names).
TARGETS = [
    {
        "pdb_id": "6BQG",
        "target_name": "5-HT2C",
        # EGFR domaine kinase + erlotinib — Stamos et al., J. Biol. Chem. 2002.
        # Référence historique pour les inhibiteurs Type I quinazoline/aminopyrimidine.
    },
    {
        "pdb_id": "7SBF",
        "target_name": "PZM21",
        # Mutant de résistance double (gatekeeper + activateur) — utile pour
        # repérer les candidats actifs spécifiquement sur la forme résistante.
    },
]

RECEPTOR_ROOT = ROOT / "data" / "receptor"


def prepare_one(pdb_id: str, target_name: str) -> None:
    target_dir = RECEPTOR_ROOT / pdb_id
    config_path = target_dir / "config.json"

    if config_path.exists():
        print(f"[{pdb_id}] déjà préparé ({config_path}) — skip.")
        return

    target_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = target_dir / f"{pdb_id}.pdb"

    print(f"[{pdb_id}] [1/4] Téléchargement depuis RCSB...")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, pdb_path)
    pdb_text = pdb_path.read_text()
    print(f"[{pdb_id}]       {len(pdb_text.splitlines())} lignes téléchargées")

    print(f"[{pdb_id}] [2/4] Calcul du centre de la poche de liaison (records SITE)...")
    center = binding_site_center(pdb_text)
    if center is None:
        print(f"[{pdb_id}]       ERREUR : aucun centre trouvé — vérifier les records SITE du PDB")
        sys.exit(1)
    print(f"[{pdb_id}]       Centre : {center}")

    print(f"[{pdb_id}] [3/4] Conversion en PDBQT (Open Babel)...")
    receptor_pdbqt = target_dir / f"{pdb_id}.pdbqt"
    ok = prepare_receptor_pdbqt(pdb_path, receptor_pdbqt)
    if not ok:
        print(f"[{pdb_id}]       ERREUR : conversion échouée — le paquet 'openbabel' (apt) est-il installé ?")
        sys.exit(1)
    print(f"[{pdb_id}]       -> {receptor_pdbqt}")

    print(f"[{pdb_id}] [4/4] Sauvegarde de la configuration...")
    config = {
        "pdb_id": pdb_id,
        "target_name": target_name,
        "center": list(center),
        "box_size": [24.0, 20.0, 24.0],
        "receptor_pdbqt": str(receptor_pdbqt.relative_to(ROOT)),
    }
    config_path.write_text(json.dumps(config, indent=2))
    print(f"[{pdb_id}]       -> {config_path}")


def main() -> None:
    for target in TARGETS:
        prepare_one(target["pdb_id"], target["target_name"])
    print(f"\nTerminé ({len(TARGETS)} cible(s) configurée(s)). "
          f"Committez data/receptor/ dans le dépôt pour que le pipeline l'utilise.")


if __name__ == "__main__":
    main()
