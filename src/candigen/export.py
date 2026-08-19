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
