#!/usr/bin/env python3
"""
Migration ponctuelle : data/hall_of_fame.json a été trouvé dans un ancien
format de schéma (clés "drug_like"/"structurally_clean", pas de "id"/"mw"/
"canonical_smiles"...), incompatible avec le MoleculeRecord actuel —
`hall_of_fame.load_hall_of_fame()` plantait dessus (TypeError).

On ne dispose de fiable, dans l'ancien format, que le champ "smiles" : ce
script recalcule tous les descripteurs à partir de là (properties +
filters), et perd donc sciemment first_seen (date de découverte
d'origine, non récupérable) et recipe (non stocké non plus à l'époque) —
mis à None plutôt que devinés. Idempotent : si le fichier est déjà au
format actuel (contient déjà "id"/"mw" etc.), ne fait rien.

Usage :
    python scripts/migrate_legacy_hall_of_fame.py --path data/hall_of_fame.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from candigen.properties import compute_descriptors, MoleculeRecord  # noqa: E402
from candigen.filters import enrich_and_filter  # noqa: E402
from candigen.evolve import fitness as compute_fitness  # noqa: E402
from candigen.hall_of_fame import save_hall_of_fame  # noqa: E402
from rdkit import Chem  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", type=Path, default=Path("data/hall_of_fame.json"))
    args = parser.parse_args()

    raw = json.loads(args.path.read_text())
    if not raw:
        print("Fichier vide, rien à migrer.")
        return
    if "id" in raw[0] and "mw" in raw[0]:
        print("Déjà au format actuel (id/mw présents), rien à faire.")
        return

    print(f"Migration de {len(raw)} entrées (ancien schéma : {sorted(raw[0].keys())})")

    records: list[MoleculeRecord] = []
    for i, entry in enumerate(raw):
        smi = entry.get("smiles")
        if not smi or Chem.MolFromSmiles(smi) is None:
            print(f"  entrée {i}: SMILES absent/invalide, ignorée")
            continue
        r = compute_descriptors(f"legacy_{i:03d}", smi)
        r.source = "legacy"
        r.notes = "Migré depuis l'ancien format de hall_of_fame.json (descripteurs recalculés, first_seen inconnu)"
        records.append(r)

    records = enrich_and_filter(records)
    for r in records:
        r.fitness = compute_fitness(r)

    save_hall_of_fame(records, args.path)
    print(f"{len(records)} entrées migrées -> {args.path}")


if __name__ == "__main__":
    main()
