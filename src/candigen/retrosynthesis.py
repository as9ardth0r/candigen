"""
candigen.retrosynthesis
========================
Rétrosynthèse assistée par ordinateur (CASP) sur les molécules conformes
au TPP, via AiZynthFinder (AstraZeneca, licence MIT, recherche arborescente
Monte Carlo guidée par un réseau de neurones —
https://github.com/MolecularAI/aizynthfinder).

Dépendance lourde (modèle + stock de précurseurs, plusieurs centaines de
Mo) — volontairement PAS dans requirements.txt (fichier séparé,
requirements-retrosynthesis.txt). Automatisée dans le workflow CI (le
modèle est mis en cache entre les runs via actions/cache, jamais committé
dans le dépôt) — voir la section dédiée du README pour le détail.

Installation et usage en local :
    pip install -r requirements-retrosynthesis.txt
    download_public_data ./aizynthfinder_data   # fourni par le paquet aizynthfinder
    python scripts/run_retrosynthesis.py --config aizynthfinder_data/config.yml

Ce que ça sort : une SÉQUENCE DE RÉACTIONS (quelle liaison casser, quels
précurseurs, jusqu'à des molécules réputées achetables) — pas un protocole
expérimental complet (réactifs exacts, solvants, températures, rendements).
La prédiction fiable des conditions réactionnelles reste un sous-problème
de recherche ouvert, hors du périmètre de ce module.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any


def load_molecules_json(path: str | Path) -> list[dict[str, Any]]:
    """Charge data/molecules.json (sortie de scripts/run_pipeline.py / candigen.export)."""
    return json.loads(Path(path).read_text())


def select_candidates(records: list[dict[str, Any]], max_n: int = 10) -> list[dict[str, Any]]:
    """
    Sélectionne les molécules conformes au TPP (`tpp_pass=True`), triées par
    fitness décroissante, plafonné à `max_n` — même logique de priorisation
    que le docking (cf. scripts/run_pipeline.py, MAX_DOCKING), pour ne pas
    lancer une recherche MCTS coûteuse sur toute la bibliothèque.
    """
    passing = [r for r in records if r.get("tpp_pass")]
    passing.sort(key=lambda r: r.get("fitness") if r.get("fitness") is not None else -99, reverse=True)
    return passing[:max_n]


def build_finder(configfile: str | Path, stock: str = "zinc", policy: str = "uspto"):
    """
    Instancie et configure un AiZynthFinder (import différé — voir
    docstring du module : dépendance lourde et optionnelle, non chargée
    tant que cette fonction n'est pas appelée). Lève ImportError si le
    paquet `aizynthfinder` n'est pas installé.
    """
    from aizynthfinder.aizynthfinder import AiZynthFinder

    finder = AiZynthFinder(configfile=str(configfile))
    finder.stock.select(stock)
    finder.expansion_policy.select(policy)
    try:
        finder.filter_policy.select(policy)
    except (KeyError, ValueError):
        pass  # policy de filtrage optionnelle, absente de certaines configs
    return finder


def search_routes(finder, smiles: str) -> dict[str, Any]:
    """
    Lance la recherche MCTS pour une molécule et retourne les routes
    trouvées ainsi que les statistiques de recherche. Ne gère pas les
    exceptions elle-même — laissé à l'appelant (scripts/run_retrosynthesis.py)
    pour ne pas interrompre un lot entier à cause d'une molécule.
    """
    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    stats = finder.extract_statistics()
    return {
        "smiles": smiles,
        "stats": stats,
        # finder.routes.dicts : liste JSON-sérialisable, un dict par route.
        # Chaque dict EST directement l'arbre de réactions (pas de clé
        # séparée "reaction_tree" ni de score par route dans cette sortie
        # brute) : noeuds "mol" (smiles, in_stock, au plus un enfant
        # "reaction" = comment cette molécule a été fabriquée) en
        # alternance avec des noeuds "reaction" (dont les enfants sont les
        # réactifs de cette étape). Voir site/js/app.js (analyzeRoute /
        # renderRetroNode) pour un exemple de parcours de cette structure.
        "routes": finder.routes.dicts,
    }


def _walk_mol_nodes(node: dict[str, Any]):
    """Générateur interne : parcourt récursivement tous les nœuds de type
    "mol" d'un arbre de route (cible, intermédiaires, précurseurs confondus)."""
    if node.get("type") == "mol" and node.get("smiles"):
        yield node
    for child in node.get("children") or []:
        yield from _walk_mol_nodes(child)


def annotate_route_names(result: dict[str, Any], delay: float = 0.25) -> None:
    """
    Enrichit chaque nœud "mol" des routes d'un résultat de search_routes()
    avec un nom IUPAC PubChem (node["name"]) — pour afficher un nom plutôt
    qu'un SMILES brut dans le dashboard (cible, intermédiaires ET
    précurseurs). Modifie `result` EN PLACE.

    Déduplique les SMILES avant d'interroger le réseau : un même
    précurseur (ex. un réactif commun) apparaît souvent dans plusieurs
    routes du même résultat — une seule requête par SMILES unique, pas une
    par occurrence.

    Ne lève jamais d'exception réseau : un échec de lookup laisse
    simplement node["name"] absent — le dashboard retombe alors sur le
    SMILES (comportement inchangé par rapport à avant cette fonction).
    """
    from .novelty import check_pubchem, compute_inchikey, fetch_pubchem_name, UNREACHABLE

    routes = result.get("routes") or []

    unique_smiles: set[str] = set()
    for route in routes:
        for mol_node in _walk_mol_nodes(route):
            unique_smiles.add(mol_node["smiles"])

    names: dict[str, str] = {}
    for smiles in unique_smiles:
        inchikey = compute_inchikey(smiles)
        if inchikey is None:
            continue
        pubchem = check_pubchem(inchikey)
        time.sleep(delay)
        if isinstance(pubchem, dict):
            name = fetch_pubchem_name(pubchem["cid"])
            time.sleep(delay)
            if isinstance(name, str) and name != UNREACHABLE:
                names[smiles] = name

    for route in routes:
        for mol_node in _walk_mol_nodes(route):
            if mol_node["smiles"] in names:
                mol_node["name"] = names[mol_node["smiles"]]


def summarize_result(result: dict[str, Any]) -> dict[str, Any]:
    """
    Résume un résultat de search_routes() pour l'injecter dans
    data/molecules.json (badge du dashboard) — juste de quoi savoir si
    une route a été trouvée, pas l'arbre complet (qui reste dans
    site/data/retrosynthesis/<id>.json, bien plus volumineux).
    """
    n_routes = len(result.get("routes") or [])
    return {
        "retrosynthesis_route_found": n_routes > 0,
        "retrosynthesis_n_routes": n_routes,
    }


def apply_summaries(
    records: list[dict[str, Any]], summaries: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    """
    Fusionne les résumés (par id, cf. summarize_result) dans la liste de
    records chargée depuis data/molecules.json — modifie et retourne la
    même liste. Les molécules non traitées ce run (id absent de
    `summaries`) gardent leur valeur précédente (ou None si jamais évaluées).
    """
    for r in records:
        summary = summaries.get(r.get("id"))
        if summary:
            r.update(summary)
    return records
