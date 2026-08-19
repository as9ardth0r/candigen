# MolGen-EGFR

Pipeline open source de génération et de criblage de candidats-médicaments
ciblant le domaine kinase d'**EGFR** (Epidermal Growth Factor Receptor),
avec un dashboard de suivi statique déployé automatiquement sur GitHub Pages.

> ⚠️ Usage : recherche & pédagogie en chimie numérique (CADD). Les molécules
> générées sont des hypothèses *in silico* non synthétisées ni testées ;
> aucune ne doit être considérée comme un candidat-médicament validé.

## Structure du dépôt

```
molgen-egfr/
├── src/molgen/            # cœur de calcul (package Python)
│   ├── properties.py      # chargement SMILES + calcul MW/LogP/TPSA/HBD/HBA
│   ├── filters.py         # TPP, Lipinski, SA score, alertes PAINS + BRENK
│   ├── generator.py       # briques : scaffolds, R-groups, assemblage SMILES
│   ├── evolve.py          # recettes + mutation atomique (algo génétique)
│   ├── hall_of_fame.py    # persistance des meilleures molécules entre runs
│   ├── docking_prep.py    # embedding 3D (ETKDGv3 + MMFF94) + export SDF
│   ├── docking.py         # docking réel : SITE→centre poche, PDBQT, Vina
│   ├── novelty.py         # vérification PubChem/ChEMBL (InChIKey exact)
│   └── export.py          # assemblage du JSON consommé par le site
├── scripts/
│   ├── run_pipeline.py     # orchestration : bootstrap ou lot évolutif du jour
│   ├── prepare_receptor.py # prépare data/receptor/ (une fois)
│   └── fetch_vendor.sh     # (re)télécharge le SA-scorer officiel RDKit
├── vendor/                 # sascorer.py + fpscores.pkl.gz (RDKit Contrib)
├── data/
│   ├── seed_molecules.smi  # 5 candidats curés
│   ├── hall_of_fame.json   # état persistant : meilleures molécules à ce jour
│   ├── explored.json       # état persistant : SMILES déjà testés
│   ├── last_run.json       # état persistant : date du dernier lot généré
│   ├── receptor/            # structure PDB + PDBQT + centre de la poche
│   ├── molecules.json      # sortie du pipeline (curés + hall of fame)
│   └── molecules.csv
├── site/                   # site statique (GitHub Pages)
│   ├── index.html
│   ├── js/app.js
│   └── data/
│       ├── molecules.json  # index léger (généré à chaque run CI)
│       └── conformers.json # blocs SDF, top-K conformes (chargé à la demande)
├── tests/
│   ├── test_properties.py
│   ├── test_evolve.py
│   ├── test_docking.py
│   └── test_novelty.py
└── .github/workflows/deploy.yml
```

## Démarrage local

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python scripts/prepare_receptor.py # optionnel : active le docking réel (une fois)
python scripts/run_pipeline.py     # régénère data/*.json/csv + site/data/*.json
python -m pytest tests/ -q         # tests unitaires
python -m http.server 8000 -d site # prévisualiser le dashboard localement
```

## Target Product Profile (TPP) — EGFR

Bornes utilisées par `TPPProfile` (`src/molgen/filters.py`), calibrées pour
un inhibiteur ATP-compétitif administré per os :

| Critère              | Cible        | Justification |
|-----------------------|-------------|----------------|
| MW                     | 250–500 Da  | Ro5 ; poche ATP kinase de taille moyenne |
| LogP                   | 1.0–4.5     | Perméabilité cellulaire sans excès de liaison protéique |
| TPSA                   | 40–100 Å²   | Compromis perméabilité passive / solubilité |
| HBD                    | ≤ 3         | Ro5, limite l'efflux (P-gp) |
| HBA                    | ≤ 9         | Ro5 |
| Liaisons rotables      | ≤ 10        | Flexibilité conformationnelle limitée (Veber) |
| SA score (Ertl)        | ≤ 4.0       | Voie de synthèse raisonnable (échelle 1–10) |
| Violations Lipinski    | ≤ 1         | Ro5 |
| QED                    | ≥ 0.3       | Drug-likeness global |
| Alertes PAINS          | 0           | Écarte les faux positifs de criblage connus |
| Alertes BRENK          | informatif  | Groupes réactifs/toxicophores connus (Brenk et al. 2008) — calculé et affiché, pas de rejet automatique par défaut (cf. note ci-dessous) |

Le pipeline calcule chaque critère avec RDKit (`Descriptors`, `Crippen`,
`rdMolDescriptors`, `QED`, `FilterCatalog`) et le SA score officiel RDKit
Contrib (Ertl & Schuffenhauer, *J. Cheminform.* 2009).

> **Sur BRENK** : le filtre est calculé pour chaque molécule
> (`record.toxicity_alerts`) et affiché dans le dashboard, mais ne fait PAS
> échouer le TPP par défaut (`TPPProfile.forbid_toxicity_alerts = False`).
> Raison testée concrètement : BRENK signale aussi des choix de conception
> légitimes et courants — le warhead acrylamide de M2 (inhibiteur covalent,
> même chimie que l'osimertinib, un médicament approuvé) et l'alcyne
> terminal de M1 (présent dans l'erlotinib) déclenchent tous deux une
> alerte. Un rejet automatique aurait donc été scientifiquement faux ici ;
> l'info reste disponible pour une revue au cas par cas. Passez le champ à
> `True` dans `TPPProfile` pour un criblage automatique plus strict.

## Boucle de découverte évolutive

Le pipeline ne se contente pas de regénérer la même bibliothèque à chaque
run — il **découvre de nouvelles molécules à chaque exécution** et fait
persister les meilleures dans le temps (`data/hall_of_fame.json`, committé
par le workflow) :

- **Premier run ("bootstrap")** : 300 combinaisons scaffold × aniline ×
  solubilisant sont tirées au sort (sur les 7 × 14 × 13 = **1274** possibles)
  et testées — volontairement une *partie* seulement du catalogue, pas la
  totalité, pour laisser de la marge aux jours suivants.
- **Runs suivants** (déclenchés par le cron quotidien, cf. CI/CD ci-dessous),
  `molgen.evolve` combine 3 mécanismes, 10 candidats chacun :
  1. **Exploration par recette** : combinaisons scaffold/aniline/solubilisant
     jamais testées, tirées au hasard dans le catalogue restant.
  2. **Mutation de recette** : on part des meilleures molécules connues et
     on change un seul de leurs 3 composants — voisinage dans le catalogue.
  3. **Mutation atomique** : un seul changement chimique local (ajouter un
     halogène/méthyle sur un cycle aromatique disponible, en retirer un, ou
     permuter un halogène) appliqué à une molécule connue (curée ou du hall
     of fame). **C'est le mécanisme qui garantit qu'il y a toujours du
     nouveau à découvrir** : contrairement à (1) et (2), cet espace n'est
     pas un catalogue fini de 1274 combinaisons — la chimie "drug-like" est
     estimée à ~10⁶⁰ molécules. Testé concrètement : même après épuisement
     total et simulé des 1274 recettes, la mutation atomique continue de
     produire des molécules inédites, jour après jour, sans jamais retomber
     à zéro (voir `tests/test_evolve.py`).
- Les 3 mécanismes partagent un seul critère de nouveauté — `data/explored.json`
  garde la trace de tous les **SMILES canoniques** déjà testés (peu importe
  leur origine), pour ne jamais retester deux fois la même molécule.
- **Fitness** = `QED - 0.05 × SA_score` (favorise drug-likeness élevé et
  synthèse facile — heuristique simple, ajustable dans `molgen/evolve.py`).
- Les molécules conformes au TPP sont fusionnées dans le hall of fame
  (déduplication par SMILES canonique, tri par fitness, plafonné à
  `HALL_OF_FAME_MAX` = 300 — les moins bonnes sont éliminées si de
  meilleures sont découvertes). Les **5 candidats curés** restent toujours
  affichés à part, jamais soumis à ce plafond, et peuvent eux aussi servir
  de parents à la mutation atomique.
- **Un seul lot par jour civil** (`data/last_run.json`) : relancer le
  pipeline plusieurs fois le même jour (ex. `workflow_dispatch` manuel) ne
  génère rien de plus après le premier run — testé concrètement, sinon la
  mutation atomique (jamais à court d'idées) ajouterait des molécules à
  chaque appel au lieu d'être neutre.

Pour aller plus loin que la mutation atomique locale (ex. changer un cycle,
fusionner deux fragments), deux leviers : curer plus de fragments dans
`ANILINE_HINGE_SUBSTITUENTS`/`SOLUBILIZING_ARMS`/`SCAFFOLD_LIBRARY`
(`src/molgen/generator.py`), ou les extraire automatiquement par
décomposition **BRICS** d'une base existante (ZINC, ChEMBL) — RDKit fournit
`Chem.BRICS` pour ça, non branché ici mais compatible avec `generator.py`.

Pour que le pipeline reste rapide et que le site déployé reste léger même
avec un hall of fame de centaines de molécules, l'embedding 3D (le plus
coûteux) et l'export du bloc SDF dans `site/data/conformers.json` sont
réservés aux `MAX_3D_EMBEDDINGS` (150 par défaut) meilleurs candidats
conformes. **Toute** la bibliothèque (descripteurs 2D, sans SDF) reste
néanmoins exportée et parcourable dans le dashboard, qui gère recherche,
filtre par conformité/provenance, tri (fitness, date de découverte, SA
score…), pagination ("Charger plus"), un badge "nouveau" sur les molécules
découvertes le jour même, un badge BRENK (nombre d'alertes) quand
pertinent, et un badge de score de docking (⚓) pour les molécules dockées
— voir section suivante.

## Docking réel (AutoDock Vina)

Contrairement au criblage 2D (TPP/SA/PAINS/BRENK, rapide et appliqué à
toute la bibliothèque), le docking calcule un score d'affinité prédit
(kcal/mol) à partir de la géométrie 3D réelle contre une structure
cristallographique d'EGFR — un signal plus proche d'une vraie interaction
protéine-ligand.

- **Structure utilisée** : PDB **1M17** (EGFR, domaine kinase, avec
  erlotinib — Stamos et al., *J. Biol. Chem.* 2002), choisie parce que
  c'est la référence historique pour les inhibiteurs réversibles de Type I
  à cœur quinazoline/aminopyrimidine, le chimiotype majoritaire de ce
  projet. Pour cribler contre un mutant de résistance (T790M/L858R),
  changer `PDB_ID` dans `scripts/prepare_receptor.py` (ex. `"4HJO"`) et
  relancer le script.
- **Préparation (une seule fois)** : `scripts/prepare_receptor.py`
  télécharge le PDB, calcule le centre de la poche de liaison à partir des
  résidus du record `SITE` (générique — fonctionne sur n'importe quelle
  structure PDB avec un ligand co-cristallisé, pas codé en dur pour 1M17),
  et convertit le récepteur en `.pdbqt` via Open Babel. Le résultat
  (`data/receptor/`) est committé et réutilisé à chaque run — pas
  retéléchargé/reconverti systématiquement.
- **Automatique en CI** : le workflow prépare le récepteur tout seul au
  premier run si `data/receptor/config.json` est absent (le runner GitHub
  a un accès internet complet). Pour lancer ça en local : `python
  scripts/prepare_receptor.py`.
- **Dans le pipeline** : `scripts/run_pipeline.py` docke les
  `MAX_DOCKING` (10 par défaut) meilleures molécules déjà embedées en 3D
  à `exhaustiveness=4` — un compromis vitesse/qualité mesuré concrètement
  sur ce projet (~30s/molécule sur un runner 2 cœurs ; ~55s à
  exhaustiveness=8). Le score est stocké dans `record.docking_score` et
  affiché dans le dashboard (badge ⚓, triable).
- **Testé de bout en bout** (`tests/test_docking.py`) : parsing des
  records `SITE`, conversion récepteur (Open Babel) et ligand (Meeko), et
  un docking réel complet — pas seulement des mocks.

```bash
python scripts/prepare_receptor.py   # une fois, télécharge + prépare data/receptor/
python scripts/run_pipeline.py       # docke automatiquement le top-10 si le récepteur existe
```

> ⚠️ Les scores de docking sont des **prédictions rapides** (scoring
> function empirique de Vina), pas des mesures. Ils permettent de classer
> relativement des molécules entre elles, pas d'estimer une affinité
> absolue fiable — pour ça, il faudrait des méthodes plus coûteuses
> (MM-GBSA, FEP) hors du périmètre de ce pipeline.

## Vérification de nouveauté (PubChem / ChEMBL)

Une molécule "générée" n'est pas forcément une molécule "nouvelle" — elle
a pu être fabriquée/publiée par quelqu'un d'autre sans qu'on le sache.
`molgen/novelty.py` compare chaque molécule conforme au TPP aux deux
grandes bases publiques de chimie, **par correspondance exacte
d'InChIKey** (pas de similarité floue — on répond juste "cette structure
exacte existe-t-elle déjà ?") :

- **PubChem** (PUG-REST) et **ChEMBL** (REST API), tous deux interrogés.
- Trois états possibles, distingués explicitement (`record.is_novel`) :
  `True` = confirmé absente des deux bases, `False` = déjà répertoriée
  (`record.pubchem_cid` et/ou `record.chembl_id` renseignés), `None` =
  pas vérifiée (une des deux bases injoignable — **jamais interprété
  comme "nouvelle"**, un bug réel rencontré et corrigé pendant le
  développement : un réseau indisponible se traduisait à tort par
  `is_novel=True` avant que `_fetch_json` distingue "404 confirmé" de
  "requête qui a échoué"). Voir `tests/test_novelty.py`, en particulier
  `test_check_novelty_unreachable_is_indeterminate_not_true`.
- Délai de 0,25s entre appels PubChem (limite de débit recommandée : 5
  req/s), et plafonné à `MAX_NOVELTY_CHECKS` (50 par défaut) molécules
  par run, priorité aux meilleures (fitness) — sinon le bootstrap (jusqu'à
  300 candidats conformes d'un coup) dépasserait largement la minute rien
  qu'en délais de limite de débit.
- Affiché dans le dashboard : badge rose "absente de PubChem/ChEMBL" si
  confirmée nouvelle, badge orange "déjà connue" avec lien direct vers la
  fiche PubChem/ChEMBL sinon — rien si pas encore vérifiée.

> ℹ️ Comme pour le téléchargement du PDB (docking), ces appels réseau ne
> peuvent être vérifiés qu'en environnement avec accès internet complet
> (CI) — les API PUG-REST/ChEMBL ne sont pas indexées par les moteurs de
> recherche, donc invérifiables depuis un environnement de développement
> à accès réseau restreint. Le format des requêtes est documenté de façon
> cohérente par plusieurs sources indépendantes ; toute la logique de
> parsing et de gestion d'erreur est testée localement avec des réponses
> construites selon ce format.

## CI/CD

`.github/workflows/deploy.yml` s'exécute :
- à chaque push sur `main` (changement de code) ;
- **tous les jours à 03:00 UTC** (`schedule: cron`) — c'est ce qui fait
  tourner la boucle de découverte quotidienne sans intervention ;
- manuellement, via l'onglet Actions → "Run workflow" (`workflow_dispatch`).

À chaque run :
1. installe Open Babel (`apt`, conversion PDB → PDBQT du récepteur) puis
   les dépendances Python,
2. prépare le récepteur EGFR (`data/receptor/`) s'il n'existe pas encore —
   seulement au tout premier run, ensuite réutilisé,
3. relance `scripts/run_pipeline.py` (bootstrap ou lot évolutif du jour,
   criblage 2D, puis docking réel du top-10 si le récepteur est prêt),
4. **committe** `data/hall_of_fame.json`, `data/explored.json`,
   `data/last_run.json`, `data/receptor/` et les fichiers `site/data/*.json`
   régénérés — sinon tout serait reperdu au run suivant, les runners
   GitHub étant éphémères. Le message de commit contient `[skip ci]` pour
   ne pas se redéclencher lui-même en boucle.
5. déploie `site/` sur GitHub Pages via `actions/deploy-pages`.

**Activation** (une fois, dans les paramètres du dépôt GitHub) :
`Settings → Pages → Source → GitHub Actions`. Le workflow a besoin de la
permission `contents: write` (déjà réglée dans `deploy.yml`) pour committer
en retour.

> ℹ️ GitHub désactive automatiquement les workflows `schedule` après 60
> jours sans **aucun push** sur le dépôt. Comme ce workflow committe
> lui-même chaque jour, le dépôt reste actif — mais si le cron s'arrête un
> jour sans raison apparente, un `workflow_dispatch` manuel suffit à le
> relancer.

## Licence

MIT — voir `LICENSE`.
