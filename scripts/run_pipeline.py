#!/usr/bin/env python3
"""Orchestration du cycle de veille : charge les molécules seed, fait
évoluer une population, exporte les résultats pour le dashboard.

Usage :
    python scripts/run_pipeline.py [--population 40] [--generations 15]
"""
from __future__ import annotations

import argparse
from pathlib import Path

from molgen_egfr import evolve, export, hall_of_fame

REPO_ROOT = Path(__file__).resolve().parent.parent


def load_seeds(path: Path) -> list[str]:
    seeds = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smiles = line.split("\t")[0]
            seeds.append(smiles)
    return seeds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path, default=REPO_ROOT / "data" / "seed_molecules.smi")
    parser.add_argument("--population", type=int, default=40)
    parser.add_argument("--generations", type=int, default=15)
    parser.add_argument("--seed", type=int, default=None, help="graine aléatoire (reproductibilité)")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "site" / "data" / "molecules.json")
    parser.add_argument("--hall-of-fame", type=Path, default=REPO_ROOT / "data" / "hall_of_fame.json")
    args = parser.parse_args()

    seeds = load_seeds(args.seeds)
    print(f"[run_pipeline] {len(seeds)} molécules seed chargées depuis {args.seeds}")

    result = evolve.run(
        seed_smiles=seeds,
        population_size=args.population,
        generations=args.generations,
        seed=args.seed,
    )

    print(f"[run_pipeline] {len(result.history)} générations calculées")
    if result.history:
        last = result.history[-1]
        print(f"[run_pipeline] fitness finale — meilleure: {last.best_fitness} | moyenne: {last.mean_fitness}")

    export.export_molecules(result, args.out)
    print(f"[run_pipeline] résultats exportés vers {args.out}")

    ranked = hall_of_fame.merge_and_save(args.hall_of_fame, result.final_population)
    print(f"[run_pipeline] hall of fame mis à jour ({len(ranked)} candidats) -> {args.hall_of_fame}")


if __name__ == "__main__":
    main()
