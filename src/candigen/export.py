"""
candigen.export
=============
Assemble les données finales dans le schéma JSON consommé par le site web.

Pour rester rapide même avec une bibliothèque de plusieurs centaines/milliers
de molécules, la sortie est scindée en deux fichiers :
  - site/data/molecules.json  : index léger (descripteurs, PAS de bloc 3D) —
                                 chargé une fois au démarrage du dashboard.
  - site/data/conformers.json : dict {id: bloc SDF}, uniquement pour les
                                 molécules conformes au TPP — chargé à la
                                 demande (lazy), la première fois que
                                 l'utilisateur ouvre le détail d'une molécule.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from .properties import MoleculeRecord


def read_target_name(root: str | Path, default: str = "EGFR") -> str:
    """
    Lit target_name depuis data/receptor/config.json (écrit par
    scripts/prepare_receptor.py) — source de vérité unique pour le nom de
    cible affiché (dashboard, badge), sans avoir à modifier le code des
    scripts d'orchestration à chaque changement de cible. Retombe sur
    `default` si le récepteur n'est pas encore préparé (bootstrap) ou si
    le fichier est absent/corrompu.
    """
    config_path = Path(root) / "data" / "receptor" / "config.json"
    if config_path.exists():
        try:
            return json.loads(config_path.read_text()).get("target_name", default)
        except (json.JSONDecodeError, OSError):
            pass
    return default


def build_site_payload(records: list[MoleculeRecord], target: str = "EGFR") -> dict:
    """Index léger : tous les descripteurs, sans les blocs SDF (voir build_conformers_payload)."""
    return {
        "target": target,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "count": len(records),
        "molecules": [r.to_dict() for r in records],
    }


def build_conformers_payload(sdf_blocks: dict[str, str]) -> dict:
    """Fichier séparé, chargé à la demande par le dashboard (cf. site/js/app.js)."""
    return dict(sdf_blocks)


def write_json(payload: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


# Alias conservé pour compatibilité avec le nom utilisé précédemment.
write_site_payload = write_json
