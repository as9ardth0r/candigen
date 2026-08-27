"""
candigen.evolve
=============
Boucle de découverte évolutive — deux mécanismes de génération, combinés
chaque jour :

  1. **Exploration par recette** (scaffold × aniline × solubilisant) :
     tant que l'espace combinatoire fixe (cf. `generator.py`) contient des
     combinaisons jamais testées, on en tire au sort.
  2. **Mutation atomique** : un seul changement chimique local (ajouter un
     halogène/méthyle, en retirer un, permuter un halogène) appliqué à une
     molécule déjà bonne (curée ou dans le hall of fame). Contrairement à
     (1), cet espace n'est PAS un catalogue fini — chaque molécule mutée
     peut, le jour suivant, être mutée à nouveau. C'est le mécanisme qui
     garantit qu'il y a toujours quelque chose de nouveau à découvrir, même
     après épuisement complet de l'espace combinatoire par recettes.

Les deux mécanismes partagent un seul et même critère de nouveauté : le
SMILES canonique du candidat ne doit jamais avoir été testé auparavant
(`explored`, cf. candigen.hall_of_fame) — peu importe comment il a été généré.
"""

from __future__ import annotations

import hashlib
import random
from dataclasses import dataclass, asdict
from itertools import product

from rdkit import Chem
from rdkit.Chem import RWMol

from .generator import (
    SCAFFOLD_LIBRARY,
    ANILINE_HINGE_SUBSTITUENTS,
    SOLUBILIZING_ARMS,
    _attach,
)
from .properties import MoleculeRecord

# --- Génération par recette (catalogue fixe scaffold × aniline × solubilisant) ---


@dataclass(frozen=True)
class Recipe:
    """Une combinaison scaffold × aniline × solubilisant, identifiée par ses 3 noms."""

    scaffold: str
    aniline: str
    solubilizer: str

    def mol_id(self) -> str:
        return f"gen_{self.scaffold}_{self.aniline}_{self.solubilizer}"

    def to_dict(self) -> dict:
        return asdict(self)


def build_smiles(recipe: Recipe) -> str | None:
    """Construit le SMILES correspondant à une recette (mêmes briques que generator.py)."""
    scaffold_smi = SCAFFOLD_LIBRARY.get(recipe.scaffold)
    aniline_smi = ANILINE_HINGE_SUBSTITUENTS.fragments.get(recipe.aniline)
    solubilizer_smi = SOLUBILIZING_ARMS.fragments.get(recipe.solubilizer)
    if not (scaffold_smi and aniline_smi and solubilizer_smi):
        return None
    step1 = _attach(scaffold_smi, 1, aniline_smi)
    if step1 is None:
        return None
    return _attach(step1, 2, solubilizer_smi)


def full_library() -> list[Recipe]:
    """Toutes les combinaisons possibles (utilisé au bootstrap et pour l'exploration fraîche)."""
    return [
        Recipe(s, a, b)
        for s, a, b in product(SCAFFOLD_LIBRARY, ANILINE_HINGE_SUBSTITUENTS.fragments, SOLUBILIZING_ARMS.fragments)
    ]


def random_recipe(rng: random.Random) -> Recipe:
    return Recipe(
        scaffold=rng.choice(list(SCAFFOLD_LIBRARY)),
        aniline=rng.choice(list(ANILINE_HINGE_SUBSTITUENTS.fragments)),
        solubilizer=rng.choice(list(SOLUBILIZING_ARMS.fragments)),
    )


# --- Mutation atomique (espace non fini : F/Cl/Br/CH3 ajoutés, retirés,
#     permutés sur un cycle aromatique disponible) ---

_HALOGENS = {"F": 9, "Cl": 17, "Br": 35}
_ADDABLE = {"F": 9, "Cl": 17, "Br": 35, "CH3": 6}


def _add_substituent(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    candidates = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetIsAromatic() and a.GetSymbol() == "C" and a.GetTotalNumHs() > 0
    ]
    if not candidates:
        return None
    idx = rng.choice(candidates)
    atomic_num = rng.choice(list(_ADDABLE.values()))
    rw = RWMol(mol)
    new_idx = rw.AddAtom(Chem.Atom(atomic_num))
    rw.AddBond(idx, new_idx, Chem.BondType.SINGLE)
    return _safe_sanitize(rw)


def _remove_substituent(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetDegree() == 1 and atom.GetSymbol() in _HALOGENS:
            candidates.append(atom.GetIdx())
        elif atom.GetDegree() == 1 and atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 3:
            nbr = atom.GetNeighbors()[0]
            if nbr.GetIsAromatic():
                candidates.append(atom.GetIdx())
    if not candidates:
        return None
    idx = rng.choice(candidates)
    rw = RWMol(mol)
    rw.RemoveAtom(idx)
    return _safe_sanitize(rw)


def _swap_halogen(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    candidates = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _HALOGENS]
    if not candidates:
        return None
    idx = rng.choice(candidates)
    current = mol.GetAtomWithIdx(idx).GetSymbol()
    options = [s for s in _HALOGENS if s != current]
    rw = RWMol(mol)
    rw.GetAtomWithIdx(idx).SetAtomicNum(_HALOGENS[rng.choice(options)])
    return _safe_sanitize(rw)


def _safe_sanitize(rw: RWMol) -> Chem.Mol | None:
    try:
        m = rw.GetMol()
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return None


_ATOM_OPS = [_add_substituent, _remove_substituent, _swap_halogen]


def mutate_atoms(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    """Applique UN changement chimique local choisi au hasard. Essaie les 3
    opérateurs dans un ordre aléatoire si le premier n'est pas applicable
    (ex. pas de cycle aromatique disponible pour `_add_substituent`)."""
    ops = list(_ATOM_OPS)
    rng.shuffle(ops)
    for op in ops:
        result = op(mol, rng)
        if result is not None:
            return result
    return None


# --- Fitness et lot quotidien ---


def fitness(record: MoleculeRecord) -> float:
    """
    Score composite simple pour classer/sélectionner les molécules :
    favorise un QED élevé (drug-likeness) et un SA score bas (synthèse facile).
    Purement heuristique — à ajuster selon vos priorités (ex. pondérer TPSA
    si la perméabilité cellulaire est prioritaire).
    """
    qed = record.qed if record.qed is not None else 0.0
    sa = record.sa_score if record.sa_score is not None else 10.0
    return round(qed - 0.05 * sa, 4)


def generate_daily_batch(
    explored: set[str],
    elite_records: list[MoleculeRecord],
    n_fresh: int,
    n_recipe_mutants: int,
    n_atom_mutants: int,
    seed: str,
) -> list[tuple[str, str, dict | None]]:
    """
    Construit le lot de candidats à tester aujourd'hui. Retourne une liste
    de (mol_id, smiles, recipe_ou_None). `explored` est l'ensemble des
    SMILES CANONIQUES déjà testés (peu importe leur origine) — le seul
    critère de nouveauté, unifié entre recettes et mutations atomiques.

      - `n_fresh`         : combinaisons scaffold/aniline/solubilisant
                            jamais testées (0 une fois le catalogue épuisé —
                            c'est attendu, pas une erreur).
      - `n_recipe_mutants`: recettes voisines des meilleures molécules
                            connues (change un seul des 3 composants).
      - `n_atom_mutants`  : mutations atomiques de molécules connues
                            (curées ou du hall of fame) — espace non fini,
                            reste productif même quand les deux précédents
                            n'ont plus rien à offrir.

    `seed` doit être une chaîne stable pour la journée (ex. la date ISO)
    afin que le lot soit reproductible si le pipeline est relancé le même jour.
    """
    rng = random.Random(seed)
    out: list[tuple[str, str, dict | None]] = []
    seen = set(explored)

    def _try_add(mol_id: str, smi: str | None, recipe: dict | None) -> bool:
        if smi is None:
            return False
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return False
        canon = Chem.MolToSmiles(mol)
        if canon in seen:
            return False
        seen.add(canon)
        out.append((mol_id, smi, recipe))
        return True

    # 1) Exploration : combinaisons neuves tirées au hasard dans le catalogue
    total_space = len(SCAFFOLD_LIBRARY) * len(ANILINE_HINGE_SUBSTITUENTS.fragments) * len(SOLUBILIZING_ARMS.fragments)
    added, attempts = 0, 0
    while added < n_fresh and attempts < total_space * 3:
        attempts += 1
        r = random_recipe(rng)
        if _try_add(r.mol_id(), build_smiles(r), r.to_dict()):
            added += 1

    # 2) Recettes voisines des meilleures molécules connues (seulement pour
    #    celles qui ONT une recette — les mutants atomiques n'en ont pas)
    elite_recipes = [Recipe(**r.recipe) for r in elite_records if r.recipe]
    if elite_recipes:
        added, attempts = 0, 0
        while added < n_recipe_mutants and attempts < n_recipe_mutants * 10:
            attempts += 1
            parent = rng.choice(elite_recipes)
            gene = rng.choice(["scaffold", "aniline", "solubilizer"])
            if gene == "scaffold":
                options = [s for s in SCAFFOLD_LIBRARY if s != parent.scaffold]
                child = Recipe(rng.choice(options), parent.aniline, parent.solubilizer)
            elif gene == "aniline":
                options = [a for a in ANILINE_HINGE_SUBSTITUENTS.fragments if a != parent.aniline]
                child = Recipe(parent.scaffold, rng.choice(options), parent.solubilizer)
            else:
                options = [s for s in SOLUBILIZING_ARMS.fragments if s != parent.solubilizer]
                child = Recipe(parent.scaffold, parent.aniline, rng.choice(options))
            if _try_add(child.mol_id(), build_smiles(child), child.to_dict()):
                added += 1

    # 3) Mutations atomiques : espace non fini, toujours productif — c'est
    #    le mécanisme qui prend le relais quand (1) et (2) sont épuisés.
    if elite_records:
        added, attempts = 0, 0
        while added < n_atom_mutants and attempts < n_atom_mutants * 15:
            attempts += 1
            parent = rng.choice(elite_records)
            parent_mol = Chem.MolFromSmiles(parent.canonical_smiles)
            if parent_mol is None:
                continue
            child_mol = mutate_atoms(parent_mol, rng)
            if child_mol is None:
                continue
            child_smiles = Chem.MolToSmiles(child_mol)
            # Basé sur le contenu (hash du SMILES canonique), pas sur le
            # compteur de boucle `attempts` seul : deploy.yml relance
            # run_pipeline.py à chaque push sur main, pas une seule fois par
            # jour comme le cron le laisserait penser — `seed` (la date) est
            # donc souvent identique entre plusieurs exécutions du même
            # jour, et un simple `attempts` remis à zéro à chaque appel peut
            # retomber sur la même valeur d'une exécution à l'autre. Un ID
            # dérivé du contenu reste unique par molécule quel que soit le
            # nombre d'exécutions dans la journée.
            content_hash = hashlib.sha1(child_smiles.encode()).hexdigest()[:8]
            mol_id = f"mut_{seed}_{attempts:04d}_{content_hash}"
            if _try_add(mol_id, child_smiles, None):
                added += 1

    return out
