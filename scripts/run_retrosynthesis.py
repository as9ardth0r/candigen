#!/usr/bin/env python3
"""
Rétrosynthèse (AiZynthFinder) sur les meilleures molécules conformes au TPP
issues de data/molecules.json — outil OPTIONNEL, à lancer manuellement en
local. PAS branché sur le workflow CI automatique (cf. docstring de
candigen.retrosynthesis pour la raison : dépendance lourde, modèle non
committable dans le dépôt).

Prérequis (une fois) :
    pip install -r requirements-retrosynthesis.txt
    download_public_data ./aizynthfinder_data

Usage :
    python scripts/run_retrosynthesis.py --config aizynthfinder_data/config.yml
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from candigen.retrosynthesis import build_finder, load_molecules_json, search_routes, select_candidates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="config.yml d'AiZynthFinder (voir download_public_data)")
    parser.add_argument("--input", default=str(ROOT / "data" / "molecules.json"), help="fichier JSON source")
    parser.add_argument("--output", default=str(ROOT / "data" / "retrosynthesis"), help="dossier de sortie")
    parser.add_argument("--max", type=int, default=10, help="nombre max de molécules traitées (top fitness)")
    parser.add_argument("--stock", default="zinc", help="nom du stock défini dans config.yml")
    parser.add_argument("--policy", default="uspto", help="nom de la policy définie dans config.yml")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = load_molecules_json(args.input)
    candidates = select_candidates(records, max_n=args.max)
    print(f"[1/2] {len(candidates)}/{len(records)} molécule(s) conformes au TPP sélectionnées (top fitness)")

    finder = build_finder(args.config, stock=args.stock, policy=args.policy)

    n_ok = 0
    for r in candidates:
        mol_id = r["id"]
        smiles = r["canonical_smiles"]
        print(f"→ {mol_id} ({smiles})")
        try:
            result = search_routes(finder, smiles)
        except Exception as exc:
            print(f"  échec : {exc}")
            continue
        out_file = output_dir / f"{mol_id}.json"
        out_file.write_text(json.dumps(result, indent=2, ensure_ascii=False))
        print(f"  {len(result['routes'])} route(s) → {out_file}")
        n_ok += 1

    print(f"[2/2] Terminé : {n_ok}/{len(candidates)} molécule(s) traitées avec succès")


if __name__ == "__main__":
    main()
