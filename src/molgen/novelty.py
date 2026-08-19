"""
molgen.novelty
===============
Vérification de nouveauté : la molécule existe-t-elle déjà dans une base
publique (PubChem, ChEMBL) ? Comparaison EXACTE par InChIKey — pas de
recherche par similarité. L'objectif est de savoir si on a affaire à une
structure déjà connue/caractérisée, pas de chercher des analogues.

Ces appels réseau ne peuvent pas être testés depuis un environnement à
accès internet restreint (les API PUG-REST/ChEMBL ne sont pas indexées par
les moteurs de recherche, contrairement à un fichier statique comme un PDB
RCSB) — le format des requêtes est documenté de façon cohérente par
plusieurs sources indépendantes (IUPAC FAIR Chemistry Cookbook, PubChemPy,
doc officielle ChEMBL), mais l'appel réseau lui-même n'a été vérifié qu'en
environnement avec accès internet complet (CI), pas ici.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request

from rdkit import Chem

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/{inchikey}.json"


def compute_inchikey(smiles: str) -> str | None:
    """Calcule l'InChIKey (identifiant canonique, comparaison exacte) d'un SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        key = Chem.MolToInchiKey(mol)
        return key if key else None
    except Exception:
        return None


class _Unreachable(Exception):
    """Levée en interne quand l'API n'a pas pu être contactée (réseau, timeout,
    erreur serveur) — DISTINCT d'une réponse 404, qui est une réponse valide
    signifiant "pas trouvé". Sans cette distinction, un réseau indisponible
    (ex. environnement de développement restreint) se traduirait à tort par
    is_novel=True (confondu avec une absence confirmée) plutôt que "indéterminé".
    """


def _fetch_json(url: str, timeout: float = 10.0) -> dict | None:
    """GET une URL. Retourne le JSON parsé si trouvé, None si 404 (réponse
    valide : pas trouvé), ou lève _Unreachable pour toute autre erreur
    (réseau, timeout, serveur en erreur) — à distinguer d'un 404."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "molgen-egfr/1.0 (research)"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise _Unreachable(f"HTTP {e.code}") from e
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        raise _Unreachable(str(e)) from e


def _parse_pubchem_response(data: dict | None) -> dict | None:
    """Extrait {'cid', 'url'} d'une réponse PUG-REST JSON — séparé du fetch réseau pour rester testable sans internet."""
    if not data:
        return None
    cids = (data.get("IdentifierList") or {}).get("CID")
    if not cids:
        return None
    cid = cids[0]
    return {"cid": cid, "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"}


def _parse_chembl_response(data: dict | None) -> dict | None:
    """Extrait {'chembl_id', 'url'} d'une réponse ChEMBL JSON — séparé du fetch réseau pour rester testable sans internet."""
    if not data:
        return None
    chembl_id = data.get("molecule_chembl_id")
    if not chembl_id:
        return None
    return {"chembl_id": chembl_id, "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/"}


# Sentinel distinct de None : "on n'a pas pu vérifier" != "on a vérifié, absent".
UNREACHABLE = "unreachable"


def check_pubchem(inchikey: str) -> dict | None | str:
    """Retourne {'cid', 'url'} si trouvé, None si confirmé absent (404), UNREACHABLE si le réseau/l'API n'a pas répondu."""
    try:
        return _parse_pubchem_response(_fetch_json(PUBCHEM_URL.format(inchikey=inchikey)))
    except _Unreachable:
        return UNREACHABLE


def check_chembl(inchikey: str) -> dict | None | str:
    """Retourne {'chembl_id', 'url'} si trouvé, None si confirmé absent (404), UNREACHABLE si le réseau/l'API n'a pas répondu."""
    try:
        return _parse_chembl_response(_fetch_json(CHEMBL_URL.format(inchikey=inchikey)))
    except _Unreachable:
        return UNREACHABLE


def check_novelty(smiles: str, delay: float = 0.25) -> dict:
    """
    Vérifie une molécule contre PubChem ET ChEMBL par InChIKey exact.
    `delay` (secondes) respecte la limite de débit recommandée par PubChem
    (5 req/s max — https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest).

    `is_novel` vaut :
      - False  si trouvée dans PubChem ET/OU ChEMBL (déjà connue),
      - True   si les DEUX bases ont répondu et confirment l'absence,
      - None   si l'une des deux bases n'a pas pu être contactée
               (indéterminé — ne pas interpréter comme "nouvelle").
    """
    inchikey = compute_inchikey(smiles)
    if inchikey is None:
        return {"inchikey": None, "is_novel": None, "pubchem": None, "chembl": None}

    pubchem = check_pubchem(inchikey)
    time.sleep(delay)
    chembl = check_chembl(inchikey)

    unreachable = pubchem == UNREACHABLE or chembl == UNREACHABLE
    found = (isinstance(pubchem, dict)) or (isinstance(chembl, dict))

    if found:
        is_novel = False
    elif unreachable:
        is_novel = None
    else:
        is_novel = True

    return {
        "inchikey": inchikey,
        "is_novel": is_novel,
        "pubchem": pubchem if isinstance(pubchem, dict) else None,
        "chembl": chembl if isinstance(chembl, dict) else None,
    }
