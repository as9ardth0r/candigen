#!/usr/bin/env python3
"""Orchestration du cycle de veille quotidien : génère un nouveau lot de
candidats (recettes fraîches + recettes voisines des meilleures connues +
mutations atomiques), calcule leurs descripteurs, les filtre contre le
TPP, fusionne les conformes dans le hall of fame persistant, et exporte
les données consommées par le dashboard.

Conçu pour tourner une fois par jour (cf. cron dans .github/workflows/deploy.yml)
— PAS un algorithme génétique classique "population/générations" relancé de
zéro à chaque fois : le hall of fame et explored.json sont l'état persistant
entre deux exécutions (les runners GitHub Actions sont éphémères).

Usage :
    python scripts/run_pipeline.py [--n-fresh 15] [--n-recipe-mutants 15] [--n-atom-mutants 10]
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rdkit import Chem  # noqa: E402

from candigen import evolve, properties, filters, hall_of_fame, export  # noqa: E402
from candigen.docking_prep import embed_3d, mol_to_sdf_block  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_fresh_seeds(seeds_path: Path, explored: set[str]) -> list[tuple[str, str]]:
    """Charge les molécules seed non encore explorées — même critère de
    nouveauté unifié que le reste du pipeline (cf. evolve.py), pour éviter
    de recalculer/refusionner les mêmes seeds à chaque run une fois
    qu'elles font déjà partie du hall of fame."""
    fresh = []
    for mol_id, smi in properties.load_smiles(seeds_path):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        if Chem.MolToSmiles(mol) not in explored:
            fresh.append((mol_id, smi))
    return fresh


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", type=Path, default=REPO_ROOT / "data" / "seed_molecules.smi")
    parser.add_argument("--n-fresh", type=int, default=15, help="nouvelles combinaisons scaffold/aniline/solubilisant")
    parser.add_argument("--n-recipe-mutants", type=int, default=15, help="recettes voisines des meilleures connues")
    parser.add_argument("--n-atom-mutants", type=int, default=10, help="mutations atomiques de molécules connues")
    parser.add_argument("--date", type=str, default=None, help="graine du jour, défaut = date ISO du jour (reproductibilité)")
    parser.add_argument("--molecules-json", type=Path, default=REPO_ROOT / "data" / "molecules.json")
    parser.add_argument("--molecules-csv", type=Path, default=REPO_ROOT / "data" / "molecules.csv")
    parser.add_argument("--hall-of-fame", type=Path, default=REPO_ROOT / "data" / "hall_of_fame.json")
    parser.add_argument("--explored", type=Path, default=REPO_ROOT / "data" / "explored.json")
    parser.add_argument("--last-run", type=Path, default=REPO_ROOT / "data" / "last_run.json")
    parser.add_argument("--site-molecules", type=Path, default=REPO_ROOT / "site" / "data" / "molecules.json")
    parser.add_argument("--site-conformers", type=Path, default=REPO_ROOT / "site" / "data" / "conformers.json")
    args = parser.parse_args()

    today = args.date or date.today().isoformat()

    hall = hall_of_fame.load_hall_of_fame(args.hall_of_fame)
    explored = hall_of_fame.load_explored(args.explored)
    elites = hall_of_fame.elite_records(hall)
    print(f"[run_pipeline] hall of fame actuel : {len(hall)} molécule(s), {len(explored)} SMILES déjà explorés")

    seed_entries = _load_fresh_seeds(args.seeds, explored)
    print(f"[run_pipeline] {len(seed_entries)} molécule(s) seed non encore explorée(s)")

    batch = evolve.generate_daily_batch(
        explored=explored,
        elite_records=elites,
        n_fresh=args.n_fresh,
        n_recipe_mutants=args.n_recipe_mutants,
        n_atom_mutants=args.n_atom_mutants,
        seed=today,
    )
    print(f"[run_pipeline] {len(batch)} nouveau(x) candidat(s) généré(s) (lot du jour)")

    recipe_by_id = {mol_id: recipe for mol_id, _smi, recipe in batch}
    all_entries = seed_entries + [(mol_id, smi) for mol_id, smi, _recipe in batch]

    records = properties.compute_batch(all_entries)
    seed_ids = {mol_id for mol_id, _ in seed_entries}
    for r in records:
        if r.id in seed_ids:
            r.source = "seed"
        elif r.id in recipe_by_id:
            r.source = "generated"
            r.recipe = recipe_by_id[r.id]
        else:
            r.source = "generated"  # mutant atomique : pas de recette, id préfixé "mut_"

    records = filters.enrich_and_filter(records)
    for r in records:
        r.fitness = evolve.fitness(r)

    passing = [r for r in records if r.tpp_pass]
    print(f"[run_pipeline] {len(passing)}/{len(records)} conforme(s) au TPP")

    merged = hall_of_fame.merge_into_hall_of_fame(hall, passing, today)
    hall_of_fame.save_hall_of_fame(merged, args.hall_of_fame)

    explored |= {r.canonical_smiles for r in records}
    hall_of_fame.save_explored(explored, args.explored)
    print(f"[run_pipeline] hall of fame mis à jour : {len(merged)} candidat(s) -> {args.hall_of_fame}")

    properties.export_json(merged, args.molecules_json)
    properties.export_csv(merged, args.molecules_csv)

    target = export.read_target_name(REPO_ROOT)
    export.write_json(export.build_site_payload(merged, target=target), args.site_molecules)

    # Conformères 3D : uniquement pour les molécules conformes au TPP (cf.
    # export.py), et seulement celles qui n'ont pas déjà un conformère
    # sauvegardé lors d'un run précédent — évite de tout ré-embarquer à
    # chaque exécution alors que le hall of fame ne fait que grossir.
    existing_conformers: dict[str, str] = {}
    if args.site_conformers.exists():
        try:
            existing_conformers = json.loads(args.site_conformers.read_text())
        except json.JSONDecodeError:
            existing_conformers = {}

    sdf_blocks = dict(existing_conformers)
    n_new_conformers = 0
    for r in merged:
        if not r.tpp_pass or r.id in sdf_blocks:
            continue
        mol = embed_3d(r.canonical_smiles, r.id)
        if mol is not None:
            sdf_blocks[r.id] = mol_to_sdf_block(mol)
            n_new_conformers += 1
    # ne garde que les conformeres des molecules toujours dans le hall of fame actuel
    merged_ids = {r.id for r in merged}
    sdf_blocks = {mol_id: sdf for mol_id, sdf in sdf_blocks.items() if mol_id in merged_ids}
    export.write_json(export.build_conformers_payload(sdf_blocks), args.site_conformers)
    print(f"[run_pipeline] conformères : {n_new_conformers} nouveau(x), {len(sdf_blocks)} au total -> {args.site_conformers}")

    args.last_run.parent.mkdir(parents=True, exist_ok=True)
    args.last_run.write_text(json.dumps({"date": today}))

    print(f"[run_pipeline] terminé — site : {args.site_molecules}")


if __name__ == "__main__":
    main()
