#!/usr/bin/env python3
"""
Migration ponctuelle : les IDs de mutants atomiques générés avant le fix
de evolve.py (mol_id = f"mut_{seed}_{attempts:04d}", sans hash de contenu)
peuvent être partagés par plusieurs molécules distinctes — deploy.yml
relance run_pipeline.py à chaque push, pas une fois par jour, donc
"mut_2026-08-27_0014" a par exemple été attribué à 5 molécules différentes.

Ce script :
1. Renomme chaque ID de cette forme en y ajoutant un hash du SMILES
   canonique (même schéma que le fix), rendant chaque ID unique par
   contenu.
2. Renomme les fichiers site/data/retrosynthesis/<ancien_id>.json vers le
   nouvel ID correspondant, en les rattachant par SMILES (champ "smiles"
   du fichier) — pas en supposant un ordre, puisque plusieurs molécules
   partagent l'ancien nom.
3. Régénère entièrement site/data/conformers.json (contrairement à la
   rétrosynthèse, un conformère 3D est bon marché à recalculer — pas la
   peine de tenter un rattachement ambigu sur des données déjà ambiguës).
4. Régénère data/molecules.json, data/molecules.csv, site/data/molecules.json
   depuis le hall_of_fame.json corrigé.

Idempotent : un ID déjà au nouveau format (se terminant par un hash hex de
8 caractères après un underscore) n'est pas retouché.

Usage :
    python scripts/dedupe_mutant_ids.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rdkit import Chem  # noqa: E402
from candigen.docking_prep import embed_3d, mol_to_sdf_block  # noqa: E402
from candigen import properties, export  # noqa: E402
from candigen.hall_of_fame import load_hall_of_fame, save_hall_of_fame  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
OLD_STYLE_ID = re.compile(r"^mut_[\d-]+_\d{4}$")  # sans suffixe hash


def canon(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


def main() -> None:
    hof_path = REPO_ROOT / "data" / "hall_of_fame.json"
    hall = load_hall_of_fame(hof_path)

    id_rename: dict[str, list[tuple[str, str]]] = {}  # ancien_id -> [(nouveau_id, canonical_smiles), ...]
    n_renamed = 0
    for r in hall:
        if OLD_STYLE_ID.match(r.id):
            old_id = r.id
            content_hash = hashlib.sha1(r.canonical_smiles.encode()).hexdigest()[:8]
            new_id = f"{old_id}_{content_hash}"
            r.id = new_id
            id_rename.setdefault(old_id, []).append((new_id, r.canonical_smiles))
            n_renamed += 1

    print(f"{n_renamed} ID(s) renommé(s) (ajout d'un hash de contenu) sur {len(hall)} molécule(s).")
    save_hall_of_fame(hall, hof_path)

    # --- rétrosynthèse : rattacher par SMILES, pas par position ---
    retro_dir = REPO_ROOT / "site" / "data" / "retrosynthesis"
    n_retro_renamed, n_retro_ambiguous = 0, 0
    for old_id, candidates in id_rename.items():
        retro_path = retro_dir / f"{old_id}.json"
        if not retro_path.exists():
            continue
        data = json.loads(retro_path.read_text())
        target_canon = canon(data.get("smiles", ""))
        matches = [new_id for new_id, smi in candidates if smi == target_canon]
        if len(matches) == 1:
            new_path = retro_dir / f"{matches[0]}.json"
            retro_path.rename(new_path)
            print(f"  rétrosynthèse : {old_id}.json -> {matches[0]}.json")
            n_retro_renamed += 1
        else:
            print(f"  ATTENTION rétrosynthèse {old_id}.json : {len(matches)} correspondance(s) "
                  f"par SMILES parmi {len(candidates)} candidat(s) — non renommé, à vérifier à la main.")
            n_retro_ambiguous += 1
    print(f"{n_retro_renamed} fichier(s) de rétrosynthèse renommé(s), {n_retro_ambiguous} ambigu(s).")

    # --- conformères : régénérés entièrement (bon marché, évite tout rattachement ambigu) ---
    sdf_blocks = {}
    for r in hall:
        if not r.tpp_pass:
            continue
        mol = embed_3d(r.canonical_smiles, r.id)
        if mol is not None:
            sdf_blocks[r.id] = mol_to_sdf_block(mol)
    export.write_json(export.build_conformers_payload(sdf_blocks), REPO_ROOT / "site" / "data" / "conformers.json")
    print(f"conformères régénérés : {len(sdf_blocks)}")

    # --- exports dérivés, régénérés depuis le hall_of_fame corrigé ---
    properties.export_json(hall, REPO_ROOT / "data" / "molecules.json")
    properties.export_csv(hall, REPO_ROOT / "data" / "molecules.csv")
    target = export.read_target_name(REPO_ROOT)
    export.write_json(export.build_site_payload(hall, target=target), REPO_ROOT / "site" / "data" / "molecules.json")
    print("data/molecules.json, data/molecules.csv, site/data/molecules.json régénérés.")


if __name__ == "__main__":
    main()
