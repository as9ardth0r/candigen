#!/usr/bin/env python3
"""
Rétrosynthèse (AiZynthFinder) sur les meilleures molécules conformes au TPP
issues de data/molecules.json. Écrit un JSON détaillé par molécule dans
data/retrosynthesis/, ET met à jour le badge "route trouvée" dans
data/molecules.json, data/molecules.csv et site/data/molecules.json
(désactivable avec --no-dashboard-update).

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

from candigen.export import build_site_payload, write_json
from candigen.properties import MoleculeRecord, export_csv, export_json
from candigen.retrosynthesis import (
    apply_summaries,
    build_finder,
    load_molecules_json,
    search_routes,
    select_candidates,
    summarize_result,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="config.yml d'AiZynthFinder (voir download_public_data)")
    parser.add_argument("--input", default=str(ROOT / "data" / "molecules.json"), help="fichier JSON source")
    parser.add_argument("--output", default=str(ROOT / "data" / "retrosynthesis"), help="dossier de sortie")
    parser.add_argument("--max", type=int, default=10, help="nombre max de molécules traitées (top fitness)")
    parser.add_argument("--stock", default="zinc", help="nom du stock défini dans config.yml")
    parser.add_argument("--policy", default="uspto", help="nom de la policy définie dans config.yml")
    parser.add_argument(
        "--no-dashboard-update", action="store_true",
        help="ne pas ré-exporter data/molecules.json + site/data/molecules.json (juste écrire data/retrosynthesis/)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = load_molecules_json(args.input)
    candidates = select_candidates(records, max_n=args.max)
    print(f"[1/3] {len(candidates)}/{len(records)} molécule(s) conformes au TPP sélectionnées (top fitness)")

    finder = build_finder(args.config, stock=args.stock, policy=args.policy)

    summaries: dict[str, dict] = {}
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
        summaries[mol_id] = summarize_result(result)
        print(f"  {len(result['routes'])} route(s) → {out_file}")
        n_ok += 1

    print(f"[2/3] Terminé : {n_ok}/{len(candidates)} molécule(s) traitées avec succès")

    if args.no_dashboard_update:
        return

    records = apply_summaries(records, summaries)
    full_records = [MoleculeRecord(**r) for r in records]
    export_json(full_records, ROOT / "data" / "molecules.json")
    export_csv(full_records, ROOT / "data" / "molecules.csv")
    write_json(build_site_payload(full_records, target="EGFR"), ROOT / "site" / "data" / "molecules.json")
    print(f"[3/3] Badge \"route trouvée\" mis à jour sur {len(summaries)} molécule(s) "
          f"(data/molecules.json, data/molecules.csv, site/data/molecules.json)")


if __name__ == "__main__":
    main()
