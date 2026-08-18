"""
molgen.generator
================
Génération combinatoire par décoration de scaffold ("scaffold hopping" léger) :
un cœur (scaffold) porteur de points d'attache `[*:n]` est combiné avec des
bibliothèques de substituants (R-groups) pour énumérer de nouveaux analogues.

Ce n'est pas un générateur deep-learning (aucun modèle n'est fourni ici),
mais une approche rule-based reproductible et sans dépendance GPU/API externe
— adaptée à un pipeline open source léger. Elle peut être remplacée plus tard
par un générateur génératif (REINVENT, MolGPT, etc.) sans changer l'interface
en aval (`filters.py` / `properties.py` consomment simplement des SMILES).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass
class RGroupLibrary:
    """Bibliothèque de substituants nommés, chacun en SMILES avec un point d'attache `[*]`."""

    name: str
    fragments: dict[str, str]  # nom -> SMILES avec un seul [*]


def _attach(core_smiles: str, attachment_label: int, fragment_smiles: str) -> str | None:
    """
    Remplace le point d'attache `[*:label]` du scaffold par un fragment `[*]-R` :
    on relie l'atome voisin du dummy du scaffold à l'atome voisin du dummy du
    fragment par une liaison simple, puis on supprime les deux dummies.

    (Note d'implémentation : `Chem.molzip` exige que les deux dummies portent
    le même atom-map number, ce qui obligerait à ré-étiqueter chaque fragment
    des bibliothèques pour chaque position — la reconstruction manuelle via
    RWMol ci-dessous est plus simple à réutiliser sur des bibliothèques de
    R-groups partagées entre plusieurs points d'attache.)
    """
    core = Chem.MolFromSmiles(core_smiles)
    frag = Chem.MolFromSmiles(fragment_smiles)
    if core is None or frag is None:
        return None

    core_dummy_idx = next(
        (a.GetIdx() for a in core.GetAtoms() if a.GetAtomicNum() == 0 and a.GetAtomMapNum() == attachment_label),
        None,
    )
    frag_dummy_idx = next((a.GetIdx() for a in frag.GetAtoms() if a.GetAtomicNum() == 0), None)
    if core_dummy_idx is None or frag_dummy_idx is None:
        return None

    core_neighbors = core.GetAtomWithIdx(core_dummy_idx).GetNeighbors()
    frag_neighbors = frag.GetAtomWithIdx(frag_dummy_idx).GetNeighbors()
    if not core_neighbors or not frag_neighbors:
        return None
    core_anchor = core_neighbors[0].GetIdx()
    frag_anchor = frag_neighbors[0].GetIdx()

    combo = Chem.RWMol(Chem.CombineMols(core, frag))
    offset = core.GetNumAtoms()
    combo.AddBond(core_anchor, frag_anchor + offset, Chem.BondType.SINGLE)
    for idx in sorted([core_dummy_idx, frag_dummy_idx + offset], reverse=True):
        combo.RemoveAtom(idx)

    try:
        mol = combo.GetMol()
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol)


def enumerate_analogs(scaffold_smiles: str, rgroups: dict[int, RGroupLibrary], id_prefix: str = "gen") -> list[tuple[str, str]]:
    """
    Énumère toutes les combinaisons de substituants sur UN scaffold multi-points
    d'attache. `rgroups` associe un numéro de point d'attache `[*:n]` à une
    `RGroupLibrary`. Retourne une liste de (id, smiles).

    Exemple :
        scaffold = "c1ccc2ncnc(N[*:1])c2c1[*:2]"
        rgroups  = {1: aniline_lib, 2: solubilizing_lib}
    """
    labels = sorted(rgroups.keys())
    fragment_choices = [list(rgroups[label].fragments.items()) for label in labels]

    out: list[tuple[str, str]] = []
    for combo in product(*fragment_choices):
        smi = scaffold_smiles
        names = []
        ok = True
        for label, (frag_name, frag_smi) in zip(labels, combo):
            new_smi = _attach(smi, label, frag_smi)
            if new_smi is None:
                ok = False
                break
            smi = new_smi
            names.append(frag_name)
        if ok:
            mol_id = f"{id_prefix}_{'_'.join(names)}"
            out.append((mol_id, smi))
    return out


def enumerate_library(
    scaffolds: dict[str, str],
    rgroups: dict[int, RGroupLibrary],
) -> list[tuple[str, str]]:
    """
    Énumère la bibliothèque combinatoire complète : chaque scaffold × chaque
    combinaison de substituants. C'est la fonction utilisée par
    `scripts/run_pipeline.py` pour construire un grand jeu de candidats à
    cribler (au lieu de partir d'une poignée de molécules choisies à la main).

    Taille de la bibliothèque = n_scaffolds × Π(taille de chaque RGroupLibrary).
    """
    out: list[tuple[str, str]] = []
    for scaffold_name, scaffold_smi in scaffolds.items():
        analogs = enumerate_analogs(scaffold_smi, rgroups, id_prefix=f"gen_{scaffold_name}")
        out.extend(analogs)
    return out


# --- Bibliothèques de substituants, inspirées des motifs classiques
#     d'inhibiteurs de kinases ATP-compétitifs. Élargir ces dictionnaires
#     (ou en ajouter dans SCAFFOLD_LIBRARY) augmente directement la taille
#     de la bibliothèque générée — aucun autre changement de code requis.

ANILINE_HINGE_SUBSTITUENTS = RGroupLibrary(
    name="aniline_hinge",
    fragments={
        "ethynylphenyl": "[*]c1cccc(C#C)c1",
        "chlorofluorophenyl": "[*]c1ccc(F)c(Cl)c1",
        "indazolyl": "[*]c1ccc2[nH]ncc2c1",
        "bromophenyl": "[*]c1cccc(Br)c1",
        "methoxyphenyl": "[*]c1ccc(OC)cc1",
        "trifluoromethylphenyl": "[*]c1cccc(C(F)(F)F)c1",
        "fluorophenyl": "[*]c1ccc(F)cc1",
        "chlorophenyl": "[*]c1cccc(Cl)c1",
        "cyanophenyl": "[*]c1cccc(C#N)c1",
        "methylphenyl": "[*]c1cccc(C)c1",
        "thiomethylphenyl": "[*]c1cccc(SC)c1",
        "difluorophenyl": "[*]c1ccc(F)c(F)c1",
        "trifluoromethoxyphenyl": "[*]c1cccc(OC(F)(F)F)c1",
        "dimethylaminophenyl": "[*]c1cccc(N(C)C)c1",
    },
)

SOLUBILIZING_ARMS = RGroupLibrary(
    name="solubilizer",
    fragments={
        "morpholinopropoxy": "[*]OCCCN1CCOCC1",
        "methylpiperazinylethoxy": "[*]OCCN1CCN(C)CC1",
        "sulfonylpiperidinylethoxy": "[*]OCCN1CCS(=O)(=O)CC1",
        "dimethylaminoethoxy": "[*]OCCN(C)C",
        "pyrrolidinylethoxy": "[*]OCCN1CCCC1",
        "azetidinylethoxy": "[*]OCCN1CCC1",
        "fluoropiperidinylethoxy": "[*]OCCN1CCC(F)CC1",
        "morpholinoethoxy": "[*]OCCN1CCOCC1",
        "hydroxyethylpiperazinylethoxy": "[*]OCCN1CCN(CCO)CC1",
        "cyclopropylaminoethoxy": "[*]OCCNC1CC1",
        "dimethylaminopropoxy": "[*]OCCCN(C)C",
        "oxetanylmethoxy": "[*]OCC1COC1",
        "tetrahydropyranyloxy": "[*]OC1CCOCC1",
    },
)

# Cœurs à charnière (hinge-binding cores) — chacun avec un point d'attache
# [*:1] pour l'aniline et [*:2] pour le bras solubilisant. Bioisostères les
# uns des autres (motif donneur/accepteur H différent face à la charnière).
SCAFFOLD_LIBRARY: dict[str, str] = {
    "quinazoline": "c1ccc2ncnc(N[*:1])c2c1[*:2]",
    "pyridopyrimidine": "c1nc(N[*:1])c2cc([*:2])cnc2n1",
    "pyrrolopyrimidine": "c1nc(N[*:1])c2ccn([*:2])c2n1",
    "imidazopyridine": "c1cc([*:2])n2c(N[*:1])cnc2c1",
    "quinolinepyrimidyl": "c1ccc2c(N[*:1])cc([*:2])nc2c1",
    "triazolopyrimidine": "c1nc(N[*:1])c2nnn([*:2])c2n1",
    "thienopyrimidine": "c1sc2ncnc(N[*:1])c2c1[*:2]",
}

# Taille de la bibliothèque combinatoire par défaut :
# 7 scaffolds × 14 anilines × 13 solubilisants = 1274 candidats générés
# (+ 5 curés, avant dédoublonnage). Agrandir un des trois dictionnaires
# ci-dessus fait grandir ce nombre sans autre changement de code.

