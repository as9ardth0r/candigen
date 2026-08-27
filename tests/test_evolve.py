import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rdkit import Chem

from candigen.evolve import (
    Recipe, build_smiles, random_recipe, mutate_atoms, full_library,
    generate_daily_batch, fitness,
)
from candigen.properties import compute_descriptors
from candigen.filters import TPPProfile, enrich_and_filter
from candigen.hall_of_fame import merge_into_hall_of_fame, elite_records, HALL_OF_FAME_MAX


def _record_from_recipe(recipe: Recipe) -> "compute_descriptors":
    smi = build_smiles(recipe)
    rec = compute_descriptors(recipe.mol_id(), smi)
    rec.recipe = recipe.to_dict()
    return rec


def test_random_recipe_builds_valid_smiles():
    rng = random.Random("unit-test")
    for _ in range(10):
        r = random_recipe(rng)
        smi = build_smiles(r)
        assert smi is not None
        assert Chem.MolFromSmiles(smi) is not None


def test_full_library_matches_expected_size():
    from candigen.generator import SCAFFOLD_LIBRARY, ANILINE_HINGE_SUBSTITUENTS, SOLUBILIZING_ARMS

    lib = full_library()
    expected = len(SCAFFOLD_LIBRARY) * len(ANILINE_HINGE_SUBSTITUENTS.fragments) * len(SOLUBILIZING_ARMS.fragments)
    assert len(lib) == expected
    assert len({(r.scaffold, r.aniline, r.solubilizer) for r in lib}) == expected  # toutes uniques


def test_mutate_atoms_produces_valid_molecule():
    rng = random.Random("atom-mutate-test")
    mol = Chem.MolFromSmiles("C#Cc1ccc(Nc2ncnc3cc(OC)c(OCCCN4CCN(C)CC4)cc23)c(F)c1")
    successes = 0
    for _ in range(30):
        child = mutate_atoms(mol, rng)
        if child is not None:
            Chem.SanitizeMol(child)  # ne doit jamais lever
            successes += 1
    assert successes > 20  # la grande majorité des tentatives doit réussir


def test_mutate_atoms_changes_the_molecule():
    rng = random.Random("atom-mutate-changes")
    mol = Chem.MolFromSmiles("C#Cc1ccc(Nc2ncnc3cc(OC)c(OCCCN4CCN(C)CC4)cc23)c(F)c1")
    original_canon = Chem.MolToSmiles(mol)
    child = mutate_atoms(mol, rng)
    assert child is not None
    assert Chem.MolToSmiles(child) != original_canon


def test_daily_batch_never_repeats_explored_smiles():
    rng = random.Random("seed-explored")
    sample_recipes = rng.sample(full_library(), 20)
    records = [_record_from_recipe(r) for r in sample_recipes]
    explored = {r.canonical_smiles for r in records}

    batch = generate_daily_batch(explored, records, n_fresh=10, n_recipe_mutants=10, n_atom_mutants=10, seed="2026-08-19")
    for mol_id, smi, recipe in batch:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert canon not in explored
    # pas de doublon interne au lot lui-même
    canon_smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for _, s, _ in batch]
    assert len(canon_smiles) == len(set(canon_smiles))


def test_daily_batch_is_deterministic_per_seed():
    batch_a = generate_daily_batch(set(), [], n_fresh=10, n_recipe_mutants=0, n_atom_mutants=0, seed="2026-08-19")
    batch_b = generate_daily_batch(set(), [], n_fresh=10, n_recipe_mutants=0, n_atom_mutants=0, seed="2026-08-19")
    assert [s for _, s, _ in batch_a] == [s for _, s, _ in batch_b]


def test_daily_batch_keeps_producing_after_recipe_catalog_exhausted():
    """Le cœur de la garantie : même le catalogue de recettes épuisé, la
    mutation atomique doit continuer à produire du nouveau, indéfiniment."""
    lib = full_library()
    explored = {Chem.MolToSmiles(Chem.MolFromSmiles(build_smiles(r))) for r in lib}
    assert len(explored) == len(lib)  # tout le catalogue est "déjà vu"

    rng = random.Random("elites-exhausted")
    elites = [_record_from_recipe(r) for r in rng.sample(lib, 15)]
    profile = TPPProfile()
    elites = enrich_and_filter(elites, profile)
    for r in elites:
        r.fitness = fitness(r)

    for day_offset in range(5):
        day = f"2026-09-{day_offset + 1:02d}"
        batch = generate_daily_batch(explored, elites, n_fresh=10, n_recipe_mutants=10, n_atom_mutants=10, seed=day)
        assert len(batch) > 0, f"aucun nouveau candidat au jour {day} alors que la mutation atomique devrait toujours en produire"
        for _, smi, recipe in batch:
            assert recipe is None  # forcément des mutants atomiques, le catalogue de recettes est vide
            canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
            explored.add(canon)


def test_atom_mutant_ids_dont_collide_across_same_day_runs():
    """Régression : deploy.yml relance run_pipeline.py à chaque push sur
    main, pas une seule fois par jour — generate_daily_batch peut donc être
    appelée plusieurs fois avec la même seed (la date) le même jour. L'ID
    des mutants atomiques était dérivé uniquement du compteur de boucle
    `attempts`, remis à zéro à chaque appel. Entre deux exécutions du même
    jour, la composition du hall of fame élite peut changer (le premier run
    y a fusionné de nouvelles molécules) : rng.choice(elite_records) tire
    alors un parent différent à la même position d'`attempts`, produisant
    un ID identique pour deux molécules différentes (observé en
    production : deux cartes 'mut_2026-08-27_0014' avec des SMILES
    distincts sur le dashboard)."""
    lib = full_library()
    explored = {Chem.MolToSmiles(Chem.MolFromSmiles(build_smiles(r))) for r in lib}
    rng = random.Random("collision-test")
    elites = [_record_from_recipe(r) for r in rng.sample(lib, 15)]
    profile = TPPProfile()
    elites = enrich_and_filter(elites, profile)
    for r in elites:
        r.fitness = fitness(r)

    # Run 1, avec la composition "elites_a" du hall of fame a cet instant.
    batch_a = generate_daily_batch(explored, elites, n_fresh=0, n_recipe_mutants=0, n_atom_mutants=15, seed="2026-08-27")

    # Run 2, plus tard le meme jour : le hall of fame a change de
    # composition entre-temps (nouvelles molecules fusionnees par le run 1,
    # simulees ici par un ordre different de la meme liste — suffisant pour
    # faire diverger rng.choice(elite_records) a la meme position).
    elites_shifted = list(reversed(elites))
    batch_b = generate_daily_batch(explored, elites_shifted, n_fresh=0, n_recipe_mutants=0, n_atom_mutants=15, seed="2026-08-27")

    ids_a = {mol_id: smi for mol_id, smi, _ in batch_a}
    ids_b = {mol_id: smi for mol_id, smi, _ in batch_b}
    # un ID partagé entre les deux runs est acceptable SI ET SEULEMENT SI
    # c'est la même molécule dans les deux cas (rejouée par coïncidence à
    # la même position) — le vrai bug, c'est un ID partagé pointant vers
    # des SMILES différents.
    for shared_id in set(ids_a) & set(ids_b):
        assert ids_a[shared_id] == ids_b[shared_id], (
            f"ID {shared_id} pointe vers deux molécules différentes selon le run : "
            f"{ids_a[shared_id]!r} vs {ids_b[shared_id]!r}"
        )


def test_hall_of_fame_merge_deduplicates_and_caps():
    profile = TPPProfile()
    rng = random.Random("hof-test")
    candidates = [_record_from_recipe(r) for r in rng.sample(full_library(), 20)]
    candidates = enrich_and_filter(candidates, profile)
    passing = [c for c in candidates if c.tpp_pass]

    hall = merge_into_hall_of_fame([], passing, today="2026-08-18", max_size=5)
    assert len(hall) <= 5
    assert len(hall) <= len(passing)
    assert all(r.fitness is not None and r.first_seen == "2026-08-18" for r in hall)

    hall2 = merge_into_hall_of_fame(hall, passing, today="2026-08-19", max_size=5)
    canon = [r.canonical_smiles for r in hall2]
    assert len(canon) == len(set(canon))


def test_elite_records_sorted_by_fitness():
    profile = TPPProfile()
    rng = random.Random("elite-test")
    candidates = [_record_from_recipe(r) for r in rng.sample(full_library(), 10)]
    candidates = enrich_and_filter(candidates, profile)
    passing = [c for c in candidates if c.tpp_pass]
    hall = merge_into_hall_of_fame([], passing, today="2026-08-18")
    elites = elite_records(hall, top_n=3)
    assert len(elites) <= 3
    fitnesses = [e.fitness for e in elites]
    assert fitnesses == sorted(fitnesses, reverse=True)
