"""
molgen.hall_of_fame
====================
Persistance de l'état entre deux exécutions du pipeline (les runners GitHub
Actions sont éphémères : sans ça, tout serait reperdu à chaque run).

Deux fichiers dans data/ :
  - hall_of_fame.json   : les meilleures molécules trouvées à ce jour
                           (conformes au TPP), plafonnées à HALL_OF_FAME_MAX,
                           classées par fitness décroissante.
  - explored.json        : SMILES canoniques déjà testés (recettes ET
                           mutations atomiques confondues), pour ne jamais
                           retester deux fois la même molécule.
"""

from __future__ import annotations

import json
from pathlib import Path

from .evolve import fitness as compute_fitness
from .properties import MoleculeRecord

HALL_OF_FAME_MAX = 300


def load_hall_of_fame(path: str | Path) -> list[MoleculeRecord]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as f:
        raw = json.load(f)
    return [MoleculeRecord(**entry) for entry in raw]


def save_hall_of_fame(records: list[MoleculeRecord], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump([r.to_dict() for r in records], f, indent=2, ensure_ascii=False)


def load_explored(path: str | Path) -> set[str]:
    path = Path(path)
    if not path.exists():
        return set()
    with path.open(encoding="utf-8") as f:
        return set(json.load(f))


def save_explored(explored: set[str], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(sorted(explored), f, indent=2, ensure_ascii=False)


def elite_records(hall: list[MoleculeRecord], top_n: int = 30) -> list[MoleculeRecord]:
    """Les meilleures molécules du hall of fame (pour servir de parents à la mutation)."""
    ranked = sorted(hall, key=lambda r: r.fitness if r.fitness is not None else -99, reverse=True)
    return ranked[:top_n]


def merge_into_hall_of_fame(
    hall: list[MoleculeRecord],
    new_passing: list[MoleculeRecord],
    today: str,
    max_size: int = HALL_OF_FAME_MAX,
) -> list[MoleculeRecord]:
    """
    Fusionne les nouvelles molécules conformes du jour dans le hall of fame,
    déduplique par SMILES canonique (garde la version la plus ancienne =
    first_seen le plus ancien), calcule la fitness, trie, et plafonne.
    """
    by_canonical: dict[str, MoleculeRecord] = {r.canonical_smiles: r for r in hall}
    for r in new_passing:
        r.fitness = compute_fitness(r)
        if r.canonical_smiles not in by_canonical:
            r.first_seen = today
            by_canonical[r.canonical_smiles] = r
        # sinon : déjà connue, on garde l'entrée existante (first_seen d'origine)

    merged = list(by_canonical.values())
    for r in merged:
        if r.fitness is None:
            r.fitness = compute_fitness(r)
    merged.sort(key=lambda r: r.fitness, reverse=True)
    return merged[:max_size]
