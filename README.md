# CandiGen

> **Plateforme de toxicologie prédictive, de veille anti-trafic et de cartographie rétrosynthétique des Nouvelles Substances Psychoactives (NPS)**

Pipeline open source de génération *in silico*, d'analyse structurale et de prédiction de voies de synthèse appliqué à la veille toxicologique et à la lutte contre le trafic de substances psychoactives.

> ⚠️ **Avertissement légal & éthique** : Ce projet est développé à des fins exclusives de **recherche académique en toxicologie**, de **sécurité publique** et de **veille sanitaire**. Les profils d'élucidation, mécanismes de dérivation et routes rétrosynthétiques générés *in silico* sont destinés aux laboratoires d'analyse médico-légale et aux autorités de contrôle afin de faciliter l'identification de marqueurs d'usage, d'anticiper les dérivés clandestins et d'appuyer les stratégies anti-trafic.

---

## Objectifs du projet

1. **Cartographie des dérivés et analogues (SARA)** : Modéliser l'espace chimique des analogues structuraux (fentanyloïdes, cathinones de synthèse, cannabinoïdes de synthèse, nitazènes, etc.) pour anticiper les molécules émergentes.
2. **Identification toxicologique prédictive** : Calculer les propriétés physico-chimiques, la lipophilie, l'affinité réceptorielle théorique et les profils d'alerte toxicologique.
3. **Caractérisation des voies de fabrication** : Appliquer la rétrosynthèse Assistée par Ordinateur (CASP) pour identifier les précurseurs clés, les réactifs de synthèse clandestine et appuyer la traçabilité de la chaîne d'approvisionnement.
4. **Diffusion de la veille** : Générer un dashboard statique interactif déployable sur GitHub Pages pour la consultation par les équipes de recherche ou d'expertise judiciaire.

---

## Structure du dépôt

```text
candigen/
├── src/candigen/            # Cœur d'analyse et de modélisation (package Python)
│   ├── properties.py      # Chargement SMILES + calcul MW/LogP/TPSA/HBD/HBA
│   ├── filters.py         # Profils toxicologiques, alerte d'amorce, PAINS & BRENK
│   ├── generator.py       # Génération d'analogues (squelettes, substituants, R-groups)
│   ├── evolve.py          # Exploration génétique de l'espace structural (SARA)
│   ├── hall_of_fame.py    # Persistance des structures émergentes identifiées
│   ├── docking_prep.py    # Conformateurs 3D (ETKDGv3 + MMFF94) & export SDF
│   ├── docking.py         # Docking prédictif récepteur/cible (PDBQT, AutoDock Vina)
│   ├── novelty.py         # Vérification bases officielles (PubChem/ChEMBL/UNODC)
│   ├── export.py          # Génération des schémas de données du dashboard
│   └── retrosynthesis.py  # Rétrosynthèse prédictive des voies d'accès (AiZynthFinder)
├── scripts/
│   ├── run_pipeline.py       # Orchestration globale du batch de veille quotidien
│   ├── prepare_receptor.py   # Préparation de la cible réceptorielle (ex. mu-opioïde, CB1/CB2)
│   ├── run_retrosynthesis.py # Analyse rétrosynthétique automatisée des cibles d'intérêt
│   └── fetch_vendor.sh       # Récupération des dépendances RDKit Contrib
├── vendor/                 # Modules d'évaluation de la facilité synthétique (SA score)
├── requirements.txt
├── requirements-retrosynthesis.txt  # Dépendances de planification rétrosynthétique (AiZynthFinder)
├── data/
│   ├── seed_molecules.smi  # Substances de référence (molécules mères surveillées)
│   ├── hall_of_fame.json   # Registre persistant des analogues d'intérêt toxicologique
│   ├── explored.json       # Registre des structures analysées
│   ├── last_run.json       # Horodatage du dernier cycle de veille
│   ├── receptor/            # Structure PDB / PDBQT du récepteur cible
│   ├── molecules.json      # Données globales générées pour la veille
│   └── molecules.csv
├── site/                   # Dashboard de veille stratégique (GitHub Pages)
│   ├── index.html
│   ├── js/app.js
│   └── data/
│       ├── molecules.json  # Index d'analyse
│       ├── retrosynthesis/ # Arbres de décomposition rétrosynthétique
│       └── conformers.json # Coordonnées 3D / Conformateurs des analogues
├── tests/                  # Validation et tests unitaires du pipeline
└── .github/workflows/deploy.yml
