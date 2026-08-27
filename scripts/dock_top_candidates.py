#!/usr/bin/env python3
"""Redocke les meilleurs candidats du hall of fame contre les récepteurs
disponibles. Ne s'exécute que sur une short-list (jamais sur toute la
population — voir rescore.py pour l'explication).

Usage :
    python scripts/dock_top_candidates.py --top 20 --target EGFR_T790M_C797S
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from molgen_egfr import hall_of_fame, rescore
from molgen_egfr.docking import TARGETS

REPO_ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hall-of-fame", type=Path, default=REPO_ROOT / "data" / "hall_of_fame.json")
    parser.add_argument("--top", type=int, default=20, help="nombre de candidats à redocker")
    parser.add_argument("--target", choices=list(TARGETS), default="EGFR_T790M_C797S")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "data" / "docking_results.json")
    args = parser.parse_args()

    available = rescore.available_targets()
    if args.target not in available:
        print(f"[dock_top_candidates] Récepteur '{args.target}' non préparé "
              f"(voir scripts/prepare_receptors.py). Cibles disponibles : {available or 'aucune'}.")
        return

    entries = hall_of_fame.load(args.hall_of_fame)
    if not entries:
        print(f"[dock_top_candidates] {args.hall_of_fame} est vide ou absent — "
              f"lancez d'abord scripts/run_pipeline.py.")
        return

    shortlist = [e["smiles"] for e in entries[: args.top]]
    print(f"[dock_top_candidates] docking de {len(shortlist)} candidats contre {args.target}...")

    results = rescore.rescore_with_docking(shortlist, args.target)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps([r.__dict__ for r in results], indent=2))
    print(f"[dock_top_candidates] résultats écrits -> {args.out}")
    if results:
        best = results[0]
        print(f"[dock_top_candidates] meilleur : {best.smiles} "
              f"(fitness {best.fitness}, docking {best.docking_score} kcal/mol)")


if __name__ == "__main__":
    main()
