#!/usr/bin/env python3
"""
Orchestration du pipeline — boucle de découverte évolutive.

Chaque exécution teste un lot de candidats DIFFÉRENT, jamais testé
auparavant, et fait persister les meilleurs d'un run à l'autre :

  - Premier run (aucun data/hall_of_fame.json) : "bootstrap" — un
    échantillon de la bibliothèque combinatoire (scaffolds × anilines ×
    solubilisants) est testé pour amorcer un hall of fame solide, SANS
    épuiser tout le catalogue d'un coup (il en faut pour les jours suivants).
  - Runs suivants (cron quotidien) : molgen.evolve combine
      1) de l'exploration dans le catalogue de recettes restant,
      2) des recettes voisines des meilleures molécules connues,
      3) des MUTATIONS ATOMIQUES (ajout/retrait/permutation d'un halogène
         ou d'un méthyle) des meilleures molécules connues — un espace non
         fini qui reste productif même après épuisement complet du
         catalogue de recettes (1274 combinaisons au total).

Tous les candidats passent par le même criblage 2D (molgen.properties +
molgen.filters : TPP, SA score, PAINS, BRENK). Ceux qui sont conformes au
TPP sont fusionnés dans data/hall_of_fame.json (plafonné, classé par
fitness). Le site consomme : les 5 candidats curés + le hall of fame.

Usage :
    python scripts/run_pipeline.py
"""

from __future__ import annotations

import json
import random
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from molgen.properties import load_smiles, compute_batch, compute_descriptors, export_json, export_csv, InvalidSMILESError
from molgen.filters import TPPProfile, enrich_and_filter
from molgen.evolve import full_library, generate_daily_batch, build_smiles, fitness as compute_fitness
from molgen.hall_of_fame import (
    load_hall_of_fame, save_hall_of_fame, load_explored, save_explored,
    elite_records, merge_into_hall_of_fame, HALL_OF_FAME_MAX,
)
from molgen.docking_prep import embed_3d, mol_to_sdf_block
from molgen.docking import dock as run_docking, prepare_ligand_pdbqt
from molgen.export import build_site_payload, build_conformers_payload, write_json
from rdkit import Chem

# Lot quotidien (runs après le bootstrap) : n_fresh recettes jamais testées
# + n_recipe_mutants recettes voisines des meilleures connues + n_atom_mutants
# mutations atomiques (espace non fini — cf. molgen/evolve.py).
N_FRESH_PER_RUN = 10
N_RECIPE_MUTANTS_PER_RUN = 10
N_ATOM_MUTANTS_PER_RUN = 10

# Le premier run n'échantillonne qu'UNE PARTIE du catalogue de recettes
# (pas les 1274 combinaisons d'un coup) — pour laisser de la marge aux
# jours suivants avant de basculer sur la mutation atomique.
BOOTSTRAP_SAMPLE_SIZE = 300

# L'embedding 3D est réservé à un top-K, pas à toute la bibliothèque
# (coûteux, et inutile pour les molécules hors du hall of fame).
MAX_3D_EMBEDDINGS = 150

# Le docking réel (AutoDock Vina) est réservé à un top-K plus restreint —
# mesuré sur ce projet : ~30s/molécule à exhaustiveness=4 (2 cœurs). Avec
# MAX_DOCKING=10, ça reste sous les ~5 min pour cette étape. Augmenter si
# votre budget CI le permet — le coût croît environ linéairement avec les
# deux paramètres.
MAX_DOCKING = 10
DOCKING_EXHAUSTIVENESS = 4

HALL_OF_FAME_PATH = ROOT / "data" / "hall_of_fame.json"
EXPLORED_PATH = ROOT / "data" / "explored.json"
LAST_RUN_PATH = ROOT / "data" / "last_run.json"
RECEPTOR_CONFIG_PATH = ROOT / "data" / "receptor" / "config.json"


def main() -> None:
    today = date.today().isoformat()

    # 1) Les 5 candidats curés — toujours recalculés, jamais soumis au
    #    plafond du hall of fame (on ne veut jamais les perdre).
    seed_path = ROOT / "data" / "seed_molecules.smi"
    curated_records = compute_batch(load_smiles(seed_path))
    for r in curated_records:
        r.source = "curated"
    print(f"[1/9] {len(curated_records)} candidats curés chargés")

    # 2) État persistant (vide au tout premier run)
    hall = load_hall_of_fame(HALL_OF_FAME_PATH)
    explored = load_explored(EXPLORED_PATH)
    is_bootstrap = not hall and not explored
    print(f"[2/9] État chargé : {len(hall)} molécules dans le hall of fame, "
          f"{len(explored)} SMILES déjà explorés"
          f"{' — BOOTSTRAP (premier run)' if is_bootstrap else ''}")

    # 3) Lot du jour — au plus UN lot par jour civil. Sans ce garde-fou, un
    #    déclenchement manuel répété (workflow_dispatch) continuerait à
    #    ajouter des molécules à chaque appel, puisque la mutation atomique
    #    n'épuise jamais vraiment l'espace — plus la propriété "idempotent
    #    si relancé le même jour" qu'on veut garder pour un comportement
    #    prévisible.
    last_run = json.loads(LAST_RUN_PATH.read_text())["date"] if LAST_RUN_PATH.exists() else None
    already_ran_today = (last_run == today) and not is_bootstrap

    profile = TPPProfile()
    curated_records = enrich_and_filter(curated_records, profile)
    for r in curated_records:
        r.fitness = compute_fitness(r)

    if already_ran_today:
        candidates = []
        print(f"[3/9] Déjà exécuté aujourd'hui ({today}) — aucun nouveau candidat, "
              f"le prochain lot sera généré au run suivant")
    elif is_bootstrap:
        rng = random.Random("bootstrap")
        recipes = rng.sample(full_library(), BOOTSTRAP_SAMPLE_SIZE)
        candidates = [(r.mol_id(), build_smiles(r), r.to_dict()) for r in recipes]
        candidates = [(i, s, r) for i, s, r in candidates if s is not None and s not in explored]
        print(f"[3/9] {len(candidates)} nouveaux candidats à tester aujourd'hui (bootstrap)")
    else:
        # Les curés + le hall of fame peuvent tous deux servir de parents
        # aux mutations atomiques (pas seulement les molécules générées).
        elites = elite_records(curated_records + hall, top_n=30)
        candidates = generate_daily_batch(
            explored=explored,
            elite_records=elites,
            n_fresh=N_FRESH_PER_RUN,
            n_recipe_mutants=N_RECIPE_MUTANTS_PER_RUN,
            n_atom_mutants=N_ATOM_MUTANTS_PER_RUN,
            seed=today,
        )
        print(f"[3/9] {len(candidates)} nouveaux candidats à tester aujourd'hui")

    # 4) Calcul des descripteurs pour le lot du jour
    batch_records = []
    for mol_id, smi, recipe in candidates:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        canon = Chem.MolToSmiles(mol)
        explored.add(canon)
        try:
            r = compute_descriptors(mol_id, smi)
        except InvalidSMILESError:
            continue
        r.source = "generated"
        r.recipe = recipe
        batch_records.append(r)
    print(f"[4/9] {len(batch_records)} candidats valides générés et décrits")

    # 5) Criblage TPP + SA score + PAINS + BRENK sur le lot du jour
    batch_records = enrich_and_filter(batch_records, profile)
    n_pass_today = sum(r.tpp_pass for r in batch_records)
    print(f"[5/9] Filtre TPP appliqué au lot du jour : {n_pass_today}/{len(batch_records)} conformes")

    # 6) Fusion dans le hall of fame (déduplication, fitness, plafond)
    new_passing = [r for r in batch_records if r.tpp_pass]
    hall = merge_into_hall_of_fame(hall, new_passing, today=today)
    print(f"[6/9] Hall of fame mis à jour : {len(hall)} molécules (plafond {HALL_OF_FAME_MAX})")

    save_hall_of_fame(hall, HALL_OF_FAME_PATH)
    save_explored(explored, EXPLORED_PATH)
    LAST_RUN_PATH.write_text(json.dumps({"date": today}), encoding="utf-8")

    # 7) Export final : curés + hall of fame accumulé
    records = curated_records + hall
    records.sort(key=lambda r: (r.source != "curated", -(r.fitness or -99)))

    passing_for_3d = [r for r in records if r.tpp_pass][:MAX_3D_EMBEDDINGS]
    sdf_blocks = {}
    mols_3d = {}
    for r in passing_for_3d:
        mol3d = embed_3d(r.canonical_smiles, mol_id=r.id)
        if mol3d is not None:
            sdf_blocks[r.id] = mol_to_sdf_block(mol3d)
            mols_3d[r.id] = mol3d
    print(f"[7/9] Conformères 3D générés pour {len(sdf_blocks)}/{len(passing_for_3d)} molécules")

    # 8) Docking réel (AutoDock Vina) — seulement si un récepteur a été
    #    préparé (scripts/prepare_receptor.py) et seulement pour le top-K
    #    des molécules déjà embedées en 3D à l'étape précédente.
    n_docked = 0
    if RECEPTOR_CONFIG_PATH.exists():
        receptor_config = json.loads(RECEPTOR_CONFIG_PATH.read_text())
        receptor_pdbqt = ROOT / receptor_config["receptor_pdbqt"]
        center = tuple(receptor_config["center"])
        box_size = tuple(receptor_config.get("box_size", (24.0, 20.0, 24.0)))
        for r in passing_for_3d[:MAX_DOCKING]:
            mol3d = mols_3d.get(r.id)
            if mol3d is None:
                continue
            ligand_pdbqt = prepare_ligand_pdbqt(mol3d)
            if ligand_pdbqt is None:
                continue
            score = run_docking(receptor_pdbqt, ligand_pdbqt, center, box_size, exhaustiveness=DOCKING_EXHAUSTIVENESS)
            if score is not None:
                r.docking_score = score
                n_docked += 1
        print(f"[8/9] Docking (AutoDock Vina, {receptor_config['pdb_id']}) : "
              f"{n_docked}/{min(len(passing_for_3d), MAX_DOCKING)} molécules dockées")
    else:
        print("[8/9] Pas de récepteur préparé (data/receptor/config.json absent) — "
              "docking ignoré. Lancer scripts/prepare_receptor.py pour l'activer.")

    export_json(records, ROOT / "data" / "molecules.json")
    export_csv(records, ROOT / "data" / "molecules.csv")

    index_payload = build_site_payload(records, target="EGFR")
    write_json(index_payload, ROOT / "site" / "data" / "molecules.json")

    conformers_payload = build_conformers_payload(sdf_blocks)
    write_json(conformers_payload, ROOT / "site" / "data" / "conformers.json")

    print(f"[9/9] Export -> {len(records)} molécules au total "
          f"({len(curated_records)} curées + {len(hall)} dans le hall of fame), "
          f"{len(sdf_blocks)} conformères 3D")


if __name__ == "__main__":
    main()
