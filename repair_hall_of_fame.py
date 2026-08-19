#!/usr/bin/env python3
"""
Script de réparation ponctuel — à lancer UNE FOIS depuis la racine du dépôt
(là où se trouve data/hall_of_fame.json), après avoir fait :

    git checkout 0f1221f -- data/hall_of_fame.json data/molecules.json data/molecules.csv site/data/molecules.json

Ce script réinjecte les scores de docking réels (retrouvés dans
data/molecules.json, qui les avait — le bug de persistance faisait qu'ils
n'étaient jamais écrits dans hall_of_fame.json) dans data/hall_of_fame.json
fraîchement restauré.
"""
import json
from pathlib import Path

hof_path = Path("data/hall_of_fame.json")
molecules_path = Path("data/molecules.json")

hall = json.loads(hof_path.read_text())
molecules = json.loads(molecules_path.read_text())

scores_by_id = {m["id"]: m["docking_score"] for m in molecules if m.get("docking_score") is not None}
print(f"Scores de docking réels trouvés dans {molecules_path} : {len(scores_by_id)}")

n_patched = 0
for m in hall:
    if m["id"] in scores_by_id:
        m["docking_score"] = scores_by_id[m["id"]]
        n_patched += 1

hof_path.write_text(json.dumps(hall, indent=2, ensure_ascii=False))
print(f"data/hall_of_fame.json réparé : {n_patched} entrées avec un score de docking réel.")
