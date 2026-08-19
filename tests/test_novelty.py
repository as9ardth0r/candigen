import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from candigen import novelty


def test_compute_inchikey_matches_known_value():
    # gefitinib — InChIKey réel, vérifié indépendamment via PubChem/Wikidata
    gefitinib = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
    assert novelty.compute_inchikey(gefitinib) == "XGALLCVXEZPNRQ-UHFFFAOYSA-N"


def test_compute_inchikey_invalid_smiles():
    assert novelty.compute_inchikey("pas_un_smiles(((") is None


def test_parse_pubchem_response_found():
    r = novelty._parse_pubchem_response({"IdentifierList": {"CID": [123631]}})
    assert r == {"cid": 123631, "url": "https://pubchem.ncbi.nlm.nih.gov/compound/123631"}


def test_parse_pubchem_response_empty():
    assert novelty._parse_pubchem_response(None) is None
    assert novelty._parse_pubchem_response({}) is None


def test_parse_chembl_response_found():
    r = novelty._parse_chembl_response({"molecule_chembl_id": "CHEMBL939"})
    assert r["chembl_id"] == "CHEMBL939"
    assert "CHEMBL939" in r["url"]


def test_parse_chembl_response_empty():
    assert novelty._parse_chembl_response(None) is None
    assert novelty._parse_chembl_response({}) is None


def test_check_novelty_known_compound():
    """Trouvée dans PubChem -> is_novel doit être False."""
    with patch.object(novelty, "check_pubchem", return_value={"cid": 123631, "url": "x"}), \
         patch.object(novelty, "check_chembl", return_value=None):
        result = novelty.check_novelty("COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", delay=0)
    assert result["is_novel"] is False
    assert result["pubchem"]["cid"] == 123631


def test_check_novelty_confirmed_absent():
    """Les deux bases répondent 'non trouvé' -> is_novel doit être True."""
    with patch.object(novelty, "check_pubchem", return_value=None), \
         patch.object(novelty, "check_chembl", return_value=None):
        result = novelty.check_novelty("C#CC1=CC=CC=C1", delay=0)
    assert result["is_novel"] is True
    assert result["pubchem"] is None and result["chembl"] is None


def test_check_novelty_unreachable_is_indeterminate_not_true():
    """Le cas qui a révélé un vrai bug : réseau injoignable ne doit JAMAIS
    donner is_novel=True (ça serait une confirmation à tort d'absence)."""
    with patch.object(novelty, "check_pubchem", return_value=novelty.UNREACHABLE), \
         patch.object(novelty, "check_chembl", return_value=novelty.UNREACHABLE):
        result = novelty.check_novelty("COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", delay=0)
    assert result["is_novel"] is None, "réseau injoignable doit donner un résultat indéterminé, pas True"


def test_check_novelty_partial_unreachable_still_indeterminate():
    """Une seule des deux bases injoignable (l'autre confirme l'absence) ->
    toujours indéterminé, pas 'nouveau' : on n'a pas la confirmation complète."""
    with patch.object(novelty, "check_pubchem", return_value=None), \
         patch.object(novelty, "check_chembl", return_value=novelty.UNREACHABLE):
        result = novelty.check_novelty("C#CC1=CC=CC=C1", delay=0)
    assert result["is_novel"] is None
