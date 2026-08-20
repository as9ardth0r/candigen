// CandiGen — logique du dashboard statique
// Aucune étape de build requise : ce fichier tourne tel quel sur GitHub Pages.

const DATA_URL = "data/molecules.json";
const CONFORMERS_URL = "data/conformers.json";

const TPP = {
  "MW (Da)": "250 – 500",
  "LogP": "1.0 – 4.5",
  "TPSA (Å²)": "40 – 100",
  "HBD max": "3",
  "HBA max": "9",
  "Rotatable bonds max": "10",
  "SA score max": "4.0 (facile→modéré)",
  "Violations Lipinski max": "1",
  "QED min": "0.3",
  "Alertes PAINS": "0 tolérée",
};

let state = { molecules: [], filtered: [], rdkitReady: null, page: 0, conformers: null, retroCache: {} };
const PAGE_SIZE = 24;

function renderTPP() {
  const list = document.getElementById("tpp-list");
  list.innerHTML = Object.entries(TPP)
    .map(([k, v]) => `<li class="flex justify-between border-b py-1" style="border-color:var(--line)"><span style="color:var(--ink-dim)">${k}</span><span class="font-mono tabular">${v}</span></li>`)
    .join("");
}

function summaryCard(label, value, accent = "var(--cyan)") {
  return `<div class="panel readout" style="--accent:${accent}">
    <p class="eyebrow">${label}</p>
    <p class="value tabular mt-1">${value}</p>
  </div>`;
}

function renderSummary(mols) {
  const nPass = mols.filter((m) => m.tpp_pass).length;
  const avgSA = (mols.reduce((s, m) => s + (m.sa_score || 0), 0) / mols.length).toFixed(2);
  const today = new Date().toISOString().slice(0, 10);
  const newToday = mols.filter((m) => m.first_seen === today).length;
  const bestFitness = mols.reduce((max, m) => (m.fitness !== null && m.fitness > max ? m.fitness : max), -Infinity);
  document.getElementById("summary-cards").innerHTML =
    summaryCard("Molécules générées", mols.length) +
    summaryCard("Conformes au TPP", `${nPass}/${mols.length}`, nPass === mols.length ? "var(--verdant)" : "var(--amber)") +
    summaryCard("Découvertes aujourd'hui", newToday, newToday > 0 ? "var(--amber)" : "var(--cyan)") +
    summaryCard("Meilleure fitness", bestFitness > -Infinity ? bestFitness.toFixed(3) : "—");
}

function renderChart(mols) {
  const ctx = document.getElementById("scatter-chart");
  const pass = mols.filter((m) => m.tpp_pass);
  const fail = mols.filter((m) => !m.tpp_pass);
  new Chart(ctx, {
    type: "scatter",
    data: {
      datasets: [
        {
          label: "Conforme TPP",
          data: pass.map((m) => ({ x: m.logp, y: m.tpsa, id: m.id })),
          backgroundColor: "#4ADE80",
        },
        {
          label: "Non conforme",
          data: fail.map((m) => ({ x: m.logp, y: m.tpsa, id: m.id })),
          backgroundColor: "#FF6B6B",
        },
      ],
    },
    options: {
      scales: {
        x: { title: { display: true, text: "LogP", color: "#82869C" }, ticks: { color: "#82869C" }, grid: { color: "#262A38" } },
        y: { title: { display: true, text: "TPSA (Å²)", color: "#82869C" }, ticks: { color: "#82869C" }, grid: { color: "#262A38" } },
      },
      plugins: {
        legend: { labels: { color: "#EDEEF3", font: { family: "'IBM Plex Mono', monospace", size: 11 } } },
        tooltip: { callbacks: { label: (c) => `${c.raw.id}: LogP ${c.raw.x}, TPSA ${c.raw.y}` } },
      },
    },
  });
}

function badge(passed) {
  return passed
    ? `<span class="tag tag--verdant border">TPP ✓</span>`
    : `<span class="tag tag--amber border">TPP ✗</span>`;
}

function sourceBadge(source) {
  return source === "curated"
    ? `<span class="tag tag--cyan border">curée</span>`
    : `<span class="tag tag--neutral border">générée</span>`;
}

function newBadge(firstSeen) {
  const today = new Date().toISOString().slice(0, 10);
  if (firstSeen !== today) return "";
  return `<span class="tag tag--amber border">nouveau</span>`;
}

function toxicityBadge(alerts) {
  if (!alerts || alerts.length === 0) return "";
  return `<span title="${alerts.join(', ')}" class="tag tag--rose border">⚠ BRENK ×${alerts.length}</span>`;
}

function dockingBadge(score) {
  if (score === null || score === undefined) return "";
  return `<span class="tag tag--cyan border">⚓ ${score} kcal/mol</span>`;
}

function retrosynthesisBadge(m) {
  if (m.retrosynthesis_route_found === true) {
    const n = m.retrosynthesis_n_routes ?? "?";
    return `<span title="${n} route(s) trouvée(s) — AiZynthFinder" class="tag tag--violet border">🧪 route trouvée</span>`;
  }
  if (m.retrosynthesis_route_found === false) {
    return `<span title="AiZynthFinder n'a trouvé aucune route vers des précurseurs achetables" class="tag tag--neutral border">🧪 aucune route</span>`;
  }
  return ""; // pas encore évaluée — on n'affiche rien plutôt que d'induire en erreur
}

function noveltyBadge(m) {
  if (m.is_novel === true) {
    return `<span class="tag tag--verdant border">✓ absente de PubChem/ChEMBL</span>`;
  }
  if (m.is_novel === false) {
    const label = m.pubchem_cid ? `PubChem CID ${m.pubchem_cid}` : (m.chembl_id || "connue");
    return `<span class="tag tag--amber border" title="Déjà répertoriée">⚠ déjà connue (${label})</span>`;
  }
  return ""; // is_novel === null : pas encore vérifié, on n'affiche rien plutôt que d'induire en erreur
}

function moleculeCard(m) {
  const title = m.chemical_name
    ? `<h3 class="text-sm font-medium leading-snug" style="color:var(--ink)" title="${m.chemical_name}">${m.chemical_name}</h3>
       <p class="text-[10px] font-mono tabular truncate mt-0.5" style="color:var(--ink-dim)">${m.id}</p>`
    : `<h3 class="font-mono text-sm tabular truncate" style="color:var(--cyan)">${m.id}</h3>`;
  return `<article data-id="${m.id}" tabindex="0" class="mol-card cursor-pointer panel p-4 transition">
    <div class="flex items-start justify-between gap-2">
      <div class="min-w-0">${title}</div>
      <div class="flex items-center gap-1 shrink-0">${newBadge(m.first_seen)}${sourceBadge(m.source)}${badge(m.tpp_pass)}</div>
    </div>
    <div class="mol-thumb h-32 my-3 flex items-center justify-center overflow-hidden" data-smiles="${m.canonical_smiles}"></div>
    <dl class="grid grid-cols-3 gap-x-2 gap-y-1 text-xs" style="color:var(--ink-dim)">
      <div>MW <span class="tabular" style="color:var(--ink)">${m.mw}</span></div>
      <div>LogP <span class="tabular" style="color:var(--ink)">${m.logp}</span></div>
      <div>TPSA <span class="tabular" style="color:var(--ink)">${m.tpsa}</span></div>
      <div>HBD <span class="tabular" style="color:var(--ink)">${m.hbd}</span></div>
      <div>HBA <span class="tabular" style="color:var(--ink)">${m.hba}</span></div>
      <div>SA <span class="tabular" style="color:var(--ink)">${m.sa_score}</span></div>
    </dl>
    <div class="mt-2 flex flex-wrap gap-1">${toxicityBadge(m.toxicity_alerts)}${dockingBadge(m.docking_score)}${noveltyBadge(m)}${retrosynthesisBadge(m)}</div>
  </article>`;
}

function renderGrid(mols) {
  state.page = 0;
  document.getElementById("molecule-grid").innerHTML = "";
  appendNextPage();
}

function appendNextPage() {
  const start = state.page * PAGE_SIZE;
  const batch = state.filtered.slice(start, start + PAGE_SIZE);
  const grid = document.getElementById("molecule-grid");
  grid.insertAdjacentHTML("beforeend", batch.map(moleculeCard).join(""));

  const newCards = Array.from(grid.querySelectorAll(".mol-card")).slice(start);
  newCards.forEach((el) => {
    drawThumbnail2D(el.querySelector(".mol-thumb"));
    el.addEventListener("click", () => openModal(el.dataset.id));
  });

  state.page += 1;
  updateLoadMoreButton();
  document.getElementById("result-count").textContent =
    `${Math.min(state.page * PAGE_SIZE, state.filtered.length)} / ${state.filtered.length} affichées`;
}

function updateLoadMoreButton() {
  const btn = document.getElementById("load-more");
  const shown = state.page * PAGE_SIZE;
  btn.classList.toggle("hidden", shown >= state.filtered.length);
}

document.getElementById("load-more").addEventListener("click", appendNextPage);

// --- RDKit.js : rendu 2D depuis SMILES (chargé une seule fois) ---
function getRDKit() {
  if (!state.rdkitReady) {
    state.rdkitReady = window.initRDKitModule();
  }
  return state.rdkitReady;
}

async function drawThumbnail2D(el, w = 200, h = 120) {
  try {
    const RDKit = await getRDKit();
    const mol = RDKit.get_mol(el.dataset.smiles);
    el.innerHTML = mol.get_svg(w, h);
    mol.delete();
  } catch (e) {
    el.innerHTML = `<span class="text-xs" style="color:var(--ink-dim)">rendu 2D indisponible</span>`;
  }
}

async function drawModal2D(smiles) {
  const el = document.getElementById("modal-2d");
  try {
    const RDKit = await getRDKit();
    const mol = RDKit.get_mol(smiles);
    el.innerHTML = mol.get_svg(320, 260);
    mol.delete();
  } catch (e) {
    el.innerHTML = `<span class="text-xs" style="color:var(--ink-dim)">rendu 2D indisponible</span>`;
  }
}

// --- 3Dmol.js : rendu du conformère 3D depuis le bloc SDF précalculé ---
// Un seul viewer (donc un seul contexte WebGL) est créé pour toute la durée
// de vie de la page, puis réutilisé (viewer.clear() + nouveau modèle) à
// chaque molécule. Créer un $3Dmol.createViewer() à chaque ouverture de
// modale, sans jamais le libérer, épuise la limite de contextes WebGL
// simultanés du navigateur (souvent ~16) après quelques dizaines d'ouvertures
// — silencieusement, sans erreur, ce qui explique un rendu qui s'arrête de
// fonctionner après une longue session sans que rien ne semble cassé.
let _viewer3d = null;

function draw3D(sdf) {
  const container = document.getElementById("modal-3d");
  if (!sdf) {
    container.innerHTML = `<div class="h-full flex items-center justify-center text-xs" style="color:var(--ink-dim)">Pas de conformère 3D</div>`;
    _viewer3d = null; // le canvas précédent vient d'être effacé par innerHTML
    return;
  }
  if (typeof $3Dmol === "undefined") {
    container.innerHTML = `<div class="h-full flex items-center justify-center text-xs text-center px-4" style="color:var(--ink-dim)">Bibliothèque 3Dmol.js non chargée (bloquée par un bloqueur de pub/VPN ?)</div>`;
    _viewer3d = null;
    return;
  }
  requestAnimationFrame(() => {
    try {
      if (!_viewer3d) {
        container.innerHTML = "";
        _viewer3d = $3Dmol.createViewer(container, { backgroundColor: "#020617" });
      } else {
        _viewer3d.clear();
      }
      const model = _viewer3d.addModel(sdf, "sdf");
      const nAtoms = model.selectedAtoms({}).length;
      if (nAtoms === 0) {
        console.error("3Dmol a chargé 0 atome depuis ce bloc SDF :", sdf);
        container.innerHTML = `<div class="h-full flex items-center justify-center text-xs text-center px-4" style="color:var(--ink-dim)">0 atome chargé — format SDF invalide (voir la console)</div>`;
        _viewer3d = null;
        return;
      }
      _viewer3d.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });
      _viewer3d.zoomTo();
      _viewer3d.render();
    } catch (e) {
      console.error("3Dmol render error:", e);
      container.innerHTML = `<div class="h-full flex items-center justify-center text-xs text-center px-4" style="color:var(--ink-dim)">Erreur de rendu 3D — voir la console (F12)</div>`;
      _viewer3d = null;
    }
  });
}

// --- conformers.json : chargé une seule fois, à la demande (lazy) ---
async function getConformers() {
  if (state.conformers === null) {
    try {
      const res = await fetch(CONFORMERS_URL);
      state.conformers = await res.json();
    } catch (e) {
      state.conformers = {};
    }
  }
  return state.conformers;
}

// --- site/data/retrosynthesis/<id>.json : un fichier par molécule, chargé
// à la demande (seulement si retrosynthesis_route_found=true, pour ne pas
// faire un fetch 404 inutile sur les molécules jamais évaluées) et mis en
// cache. Chemin relatif à site/index.html — résout vers site/data/retrosynthesis/
// une fois déployé sur Pages (qui ne sert que le contenu de site/).
const RETRO_URL_BASE = "data/retrosynthesis/";

async function getRetroResult(id) {
  if (state.retroCache[id] !== undefined) return state.retroCache[id];
  try {
    const res = await fetch(`${RETRO_URL_BASE}${id}.json`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    state.retroCache[id] = await res.json();
  } catch (e) {
    state.retroCache[id] = null;
  }
  return state.retroCache[id];
}

// Chaque route (result.routes[i]) est elle-même l'arbre de réactions, racine
// = molécule cible (pas de champ "score" par route dans la sortie brute
// d'AiZynthFinder — cf. candigen.retrosynthesis) : type "mol" (avec au plus
// un enfant "reaction" = comment elle a été fabriquée) en alternance avec
// type "reaction" (dont les enfants sont les réactifs de cette étape).
function analyzeRoute(root) {
  let steps = 0, precursors = 0, inStock = 0;
  (function walk(node) {
    if (node.type === "reaction") steps += 1;
    if (node.type === "mol") {
      const isLeaf = !(node.children || []).some((c) => c.type === "reaction");
      if (isLeaf) {
        precursors += 1;
        if (node.in_stock) inStock += 1;
      }
    }
    (node.children || []).forEach(walk);
  })(root);
  return { steps, precursors, inStock };
}

// Pas de champ "score" par route dans la sortie brute d'AiZynthFinder (cf.
// candigen.retrosynthesis) — on recommande donc la route la plus PRATIQUE
// à exécuter : le plus de précurseurs déjà en stock, et à égalité, le
// moins d'étapes. Purement indicatif (affiché avec ⭐), pas une garantie
// de faisabilité chimique supérieure aux autres routes proposées.
function pickBestRouteIndex(analyses) {
  let best = 0;
  for (let i = 1; i < analyses.length; i++) {
    const a = analyses[i], b = analyses[best];
    if (a.inStock > b.inStock || (a.inStock === b.inStock && a.steps < b.steps)) {
      best = i;
    }
  }
  return best;
}

function renderRetroNode(node) {
  if (node.type === "reaction") {
    const reactants = (node.children || []).map(renderRetroNode).join("");
    return `<div class="pl-4 ml-3 mt-1 retro-branch">
      <p class="text-[10px] mb-1" style="color:var(--ink-dim)">↓ réaction</p>
      ${reactants}
    </div>`;
  }
  // node.type === "mol"
  const reactionChild = (node.children || []).find((c) => c.type === "reaction");
  const stock = node.in_stock
    ? `<span class="tag tag--verdant border shrink-0">en stock</span>`
    : `<span class="tag tag--amber border shrink-0">à synthétiser</span>`;
  const label = node.name
    ? `<div class="text-[11px] leading-snug min-w-0 truncate" style="color:var(--ink)" title="${node.smiles}">${node.name}</div>`
    : `<div class="text-[10px] font-mono break-all min-w-0" style="color:var(--ink-dim)">${node.smiles}</div>`;
  return `<div class="flex items-center gap-2 mt-1">
      <div class="retro-mol-thumb h-14 w-20 shrink-0 flex items-center justify-center overflow-hidden" data-smiles="${node.smiles}"></div>
      ${label}
      ${stock}
    </div>${reactionChild ? renderRetroNode(reactionChild) : ""}`;
}

function renderRetroSection(result) {
  const select = document.getElementById("modal-retro-route-select");
  const statsEl = document.getElementById("modal-retro-stats");
  const treeEl = document.getElementById("modal-retro-tree");

  if (!result || !result.routes || result.routes.length === 0) {
    select.innerHTML = "";
    statsEl.textContent = "";
    treeEl.innerHTML = `<p class="text-xs" style="color:var(--ink-dim)">Détail de la route indisponible (site/data/retrosynthesis/ pas encore régénéré pour cette molécule).</p>`;
    return;
  }

  const analyses = result.routes.map(analyzeRoute);
  const bestIndex = pickBestRouteIndex(analyses);
  select.innerHTML = result.routes
    .map((_, i) => `<option value="${i}">${i === bestIndex ? "⭐ " : ""}Route ${i + 1} — ${analyses[i].steps} étape(s)</option>`)
    .join("");

  const showRoute = (i) => {
    const a = analyses[i];
    statsEl.textContent = `${a.steps} étape(s) · ${a.inStock}/${a.precursors} précurseur(s) déjà en stock`;
    treeEl.innerHTML = renderRetroNode(result.routes[i]);
    treeEl.querySelectorAll(".retro-mol-thumb").forEach((el) => drawThumbnail2D(el, 90, 60));
  };
  select.onchange = () => showRoute(Number(select.value));
  showRoute(bestIndex);
}

async function openModal(id) {
  const m = state.molecules.find((mol) => mol.id === id);
  if (!m) return;
  document.getElementById("modal-title").textContent = m.chemical_name || m.id;
  document.getElementById("modal-props").innerHTML = [
    ["ID", m.id],
    ["MW", m.mw], ["LogP", m.logp], ["TPSA", m.tpsa], ["HBD", m.hbd],
    ["HBA", m.hba], ["RotB", m.rotatable_bonds], ["SA score", m.sa_score], ["QED", m.qed],
    ["Fitness", m.fitness ?? "—"], ["Découverte le", m.first_seen ?? "—"],
    ["Docking (kcal/mol)", m.docking_score ?? "—"],
    ["Rétrosynthèse", m.retrosynthesis_route_found === true ? `${m.retrosynthesis_n_routes} route(s)` : (m.retrosynthesis_route_found === false ? "aucune route" : "—")],
  ].map(([k, v]) => `<div class="panel-2 rounded px-2 py-1"><span style="color:var(--ink-dim)">${k}:</span> <span class="font-mono">${v}</span></div>`).join("");
  const links = [];
  if (m.pubchem_cid) links.push(`<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${m.pubchem_cid}" target="_blank" rel="noopener" style="color:var(--cyan)" class="hover:underline">PubChem CID ${m.pubchem_cid}</a>`);
  if (m.chembl_id) links.push(`<a href="https://www.ebi.ac.uk/chembl/compound_report_card/${m.chembl_id}/" target="_blank" rel="noopener" style="color:var(--cyan)" class="hover:underline">${m.chembl_id}</a>`);
  document.getElementById("modal-notes").innerHTML =
    `${m.formula} — ${m.notes}` +
    (m.toxicity_alerts && m.toxicity_alerts.length ? ` | Alertes BRENK : ${m.toxicity_alerts.join(", ")}` : "") +
    (links.length ? ` | Déjà répertoriée : ${links.join(", ")}` : "") +
    (m.is_novel === true ? ` | Absente de PubChem/ChEMBL (vérifié)` : "");
  document.getElementById("modal").classList.remove("hidden");
  document.getElementById("modal").classList.add("flex");

  drawModal2D(m.canonical_smiles);
  document.getElementById("modal-3d").innerHTML = `<div class="h-full flex items-center justify-center text-xs" style="color:var(--ink-dim)">Chargement du conformère…</div>`;
  const conformers = await getConformers();
  draw3D(conformers[m.id]);

  const retroSection = document.getElementById("modal-retro");
  if (m.retrosynthesis_route_found === true) {
    retroSection.classList.remove("hidden");
    document.getElementById("modal-retro-tree").innerHTML = `<p class="text-xs" style="color:var(--ink-dim)">Chargement…</p>`;
    const result = await getRetroResult(m.id);
    renderRetroSection(result);
  } else {
    retroSection.classList.add("hidden");
  }
}

document.getElementById("modal-close").addEventListener("click", () => {
  document.getElementById("modal").classList.add("hidden");
  document.getElementById("modal").classList.remove("flex");
});

// --- Filtres / recherche / tri ---
function applyFilters() {
  const q = document.getElementById("search-input").value.toLowerCase();
  const passFilter = document.getElementById("filter-pass").value;
  const sourceFilter = document.getElementById("filter-source").value;
  const sortBy = document.getElementById("sort-by").value;

  let out = state.molecules.filter((m) => {
    const matchesQuery = m.id.toLowerCase().includes(q) || m.smiles.toLowerCase().includes(q) || (m.chemical_name || "").toLowerCase().includes(q);
    const matchesPass = passFilter === "all" || (passFilter === "pass" ? m.tpp_pass : !m.tpp_pass);
    const matchesSource = sourceFilter === "all" || m.source === sourceFilter;
    return matchesQuery && matchesPass && matchesSource;
  });
  const descendingFields = new Set(["fitness", "first_seen"]);
  out.sort((a, b) => {
    let va = a[sortBy], vb = b[sortBy];
    if (va == null && vb == null) return 0;
    if (va == null) return 1;  // valeurs manquantes (ex. curées sans first_seen) en dernier
    if (vb == null) return -1;
    const cmp = typeof va === "string" ? va.localeCompare(vb) : va - vb;
    return descendingFields.has(sortBy) ? -cmp : cmp;
  });
  state.filtered = out;
  renderGrid(out);
}

["search-input", "filter-pass", "filter-source", "sort-by"].forEach((id) =>
  document.getElementById(id).addEventListener("input", applyFilters)
);

// --- Chargement des données ---
async function main() {
  renderTPP();
  const res = await fetch(DATA_URL);
  const payload = await res.json();
  state.molecules = payload.molecules;
  document.getElementById("generated-at").textContent = `${payload.target} · généré le ${new Date(payload.generated_at).toLocaleString("fr-FR")}`;
  document.getElementById("page-subtitle").textContent = `Pipeline CADD open source — génération & criblage d'inhibiteurs de ${payload.target}`;
  document.title = `CandiGen (${payload.target}) — Suivi du pipeline de génération de molécules`;
  renderSummary(payload.molecules);
  renderChart(payload.molecules);
  applyFilters();
}

main().catch((err) => {
  console.error(err);
  document.getElementById("molecule-grid").innerHTML =
    `<p class="col-span-full" style="color:var(--rose)">Erreur de chargement de data/molecules.json — ${err.message}</p>`;
});
