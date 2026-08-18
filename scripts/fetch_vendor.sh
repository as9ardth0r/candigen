#!/usr/bin/env bash
# Récupère le script officiel RDKit Contrib pour le SA score (Ertl & Schuffenhauer).
# Le dépôt inclut déjà une copie dans vendor/ ; ce script sert à la remettre à jour.
set -euo pipefail
mkdir -p vendor
BASE="https://raw.githubusercontent.com/rdkit/rdkit/master/Contrib/SA_Score"
curl -sL -o vendor/sascorer.py "$BASE/sascorer.py"
curl -sL -o vendor/fpscores.pkl.gz "$BASE/fpscores.pkl.gz"
echo "vendor/sascorer.py et vendor/fpscores.pkl.gz mis à jour."
