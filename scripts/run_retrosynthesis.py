#!/usr/bin/env python3
"""
Rétrosynthèse (AiZynthFinder) sur les meilleures molécules conformes au TPP
issues de data/molecules.json. Écrit un JSON détaillé par molécule dans
site/data/retrosynthesis/ (comme site/data/conformers.json : le dashboard
GitHub Pages ne sert QUE le contenu de site/, pas data/ à la racine du
dépôt), ET met à jour le badge "route trouvée" dans
data/molecules.json, data/molecules.csv et site/data/molecules.json
(désactivable avec --no-dashboard-update).

Une molécule qui a déjà un fichier de sortie n'est PAS recalculée (le top
TPP change peu d'un jour à l'autre — sans ça on relance une recherche MCTS
coûteuse chaque jour pour un résultat identique). --force pour forcer le
recalcul de tout le lot sélectionné malgré le cache.

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

from candigen.export import build_site_payload, read_target_name, write_json
from candigen.properties import MoleculeRecord, export_csv, export_json
from candigen.retrosynthesis import (
    annotate_route_names,
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
    parser.add_argument("--output", default=str(ROOT / "site" / "data" / "retrosynthesis"), help="dossier de sortie")
    parser.add_argument("--max", type=int, default=10, help="nombre max de molécules traitées (top fitness)")
    parser.add_argument("--stock", default="zinc", help="nom du stock défini dans config.yml")
    parser.add_argument("--policy", default="uspto", help="nom de la policy définie dans config.yml")
    parser.add_argument(
        "--force", action="store_true",
        help="recalcule même les molécules qui ont déjà un fichier de sortie (ignore le cache)",
    )
    parser.add_argument(
        "--skip-names", action="store_true",
        help="ne pas récupérer les noms PubChem des précurseurs/intermédiaires (plus rapide, SMILES bruts affichés)",
    )
    parser.add_argument(
        "--no-dashboard-update", action="store_true",
        help="ne pas ré-exporter data/molecules.json + site/data/molecules.json (juste écrire site/data/retrosynthesis/)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = load_molecules_json(args.input)
    candidates = select_candidates(records, max_n=args.max)
    print(f"[1/3] {len(candidates)}/{len(records)} molécule(s) conformes au TPP sélectionnées (top fitness)")

    finder = None  # construit paresseusement — seulement si une molécule n'a pas déjà de résultat en cache
    summaries: dict[str, dict] = {}
    n_new, n_cached = 0, 0

    for r in candidates:
        mol_id = r["id"]
        smiles = r["canonical_smiles"]
        out_file = output_dir / f"{mol_id}.json"

        if out_file.exists() and not args.force:
            try:
                result = json.loads(out_file.read_text())
                summaries[mol_id] = summarize_result(result)
                n_cached += 1
                continue
            except (json.JSONDecodeError, OSError) as exc:
                print(f"  {mol_id} : fichier en cache illisible ({exc}) — recalcul")

        if finder is None:
            finder = build_finder(args.config, stock=args.stock, policy=args.policy)

        print(f"→ {mol_id} ({smiles})")
        try:
            result = search_routes(finder, smiles)
        except Exception as exc:
            print(f"  échec : {exc}")
            continue
        if not args.skip_names:
            annotate_route_names(result)
        out_file.write_text(json.dumps(result, indent=2, ensure_ascii=False))
        summaries[mol_id] = summarize_result(result)
        print(f"  {len(result['routes'])} route(s) → {out_file}")
        n_new += 1

    print(f"[2/3] Terminé : {n_new} nouvelle(s) recherche(s), {n_cached} déjà en cache "
          f"(sur {len(candidates)} molécule(s) sélectionnée(s))")

    if args.no_dashboard_update:
        return

    records = apply_summaries(records, summaries)
    full_records = [MoleculeRecord(**r) for r in records]
    export_json(full_records, ROOT / "data" / "molecules.json")
    export_csv(full_records, ROOT / "data" / "molecules.csv")
    write_json(build_site_payload(full_records, target=read_target_name(ROOT)), ROOT / "site" / "data" / "molecules.json")
    print(f"[3/3] Badge \"route trouvée\" mis à jour sur {len(summaries)} molécule(s) "
          f"(data/molecules.json, data/molecules.csv, site/data/molecules.json)")


if __name__ == "__main__":
    main()
