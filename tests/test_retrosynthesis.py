import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from candigen.retrosynthesis import apply_summaries, build_finder, select_candidates, summarize_result


def test_select_candidates_filters_and_sorts_by_fitness():
    records = [
        {"id": "a", "tpp_pass": True, "fitness": 0.5},
        {"id": "b", "tpp_pass": False, "fitness": 0.9},  # exclu : ne passe pas le TPP
        {"id": "c", "tpp_pass": True, "fitness": 0.8},
    ]
    selected = select_candidates(records, max_n=10)
    assert [r["id"] for r in selected] == ["c", "a"]


def test_select_candidates_respects_max_n():
    records = [{"id": str(i), "tpp_pass": True, "fitness": float(i)} for i in range(20)]
    selected = select_candidates(records, max_n=5)
    assert len(selected) == 5
    assert selected[0]["id"] == "19"  # meilleure fitness en premier


def test_select_candidates_handles_missing_fitness():
    records = [
        {"id": "a", "tpp_pass": True, "fitness": None},
        {"id": "b", "tpp_pass": True, "fitness": 1.0},
    ]
    selected = select_candidates(records, max_n=10)
    assert selected[0]["id"] == "b"


def test_build_finder_fails_gracefully_without_aizynthfinder():
    """Import différé : si `aizynthfinder` n'est pas installé (dépendance
    optionnelle, cf. requirements-retrosynthesis.txt), l'échec doit rester
    localisé à cet appel — pas planter le chargement du module candigen."""
    try:
        import aizynthfinder  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError):
            build_finder("nonexistent_config.yml")
    else:
        pytest.skip("aizynthfinder est installé — ce test ne vérifie que l'absence gracieuse")


def test_summarize_result_route_found():
    result = {"routes": [{"score": 0.9}, {"score": 0.7}]}
    summary = summarize_result(result)
    assert summary == {"retrosynthesis_route_found": True, "retrosynthesis_n_routes": 2}


def test_summarize_result_no_route():
    summary = summarize_result({"routes": []})
    assert summary == {"retrosynthesis_route_found": False, "retrosynthesis_n_routes": 0}


def test_apply_summaries_updates_matching_records_only():
    records = [
        {"id": "a", "fitness": 0.5},
        {"id": "b", "fitness": 0.8},
    ]
    summaries = {"a": {"retrosynthesis_route_found": True, "retrosynthesis_n_routes": 3}}
    updated = apply_summaries(records, summaries)
    assert updated[0]["retrosynthesis_route_found"] is True
    assert updated[0]["retrosynthesis_n_routes"] == 3
    assert "retrosynthesis_route_found" not in updated[1]  # non traitée ce run — inchangée
