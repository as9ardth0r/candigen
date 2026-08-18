// MolGen-EGFR — logique du dashboard statique
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

let state = { molecules: [], filtered: [], rdkitReady: null, page: 0, conformers: null };
const PAGE_SIZE = 24;

function renderTPP() {
  const list = document.getElementById("tpp-list");
  list.innerHTML = Object.entries(TPP)
    .map(([k, v]) => `<li class="flex justify-between border-b border-slate-800 py-1"><span class="text-slate-400">${k}</span><span class="font-mono">${v}</span></li>`)
    .join("");
}

function summaryCard(label, value, accent = "text-slate-100") {
  return `<div class="bg-slate-900 border border-slate-800 rounded-xl p-4">
    <p class="text-xs text-slate-400">${label}</p>
    <p class="text-2xl font-semibold ${accent}">${value}</p>
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
    summaryCard("Conformes au TPP", `${nPass}/${mols.length}`, nPass === mols.length ? "text-emerald-400" : "text-amber-400") +
    summaryCard("Découvertes aujourd'hui", newToday, newToday > 0 ? "text-emerald-400" : "text-slate-100") +
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
          backgroundColor: "#34d399",
        },
        {
          label: "Non conforme",
          data: fail.map((m) => ({ x: m.logp, y: m.tpsa, id: m.id })),
          backgroundColor: "#f87171",
        },
      ],
    },
    options: {
      scales: {
        x: { title: { display: true, text: "LogP", color: "#94a3b8" }, ticks: { color: "#94a3b8" }, grid: { color: "#1e293b" } },
        y: { title: { display: true, text: "TPSA (Å²)", color: "#94a3b8" }, ticks: { color: "#94a3b8" }, grid: { color: "#1e293b" } },
      },
      plugins: {
        legend: { labels: { color: "#cbd5e1" } },
        tooltip: { callbacks: { label: (c) => `${c.raw.id}: LogP ${c.raw.x}, TPSA ${c.raw.y}` } },
      },
    },
  });
}

function badge(passed) {
  return passed
    ? `<span class="text-xs px-2 py-0.5 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/30">TPP ✓</span>`
    : `<span class="text-xs px-2 py-0.5 rounded-full bg-amber-500/10 text-amber-400 border border-amber-500/30">TPP ✗</span>`;
}

function sourceBadge(source) {
  return source === "curated"
    ? `<span class="text-[10px] px-1.5 py-0.5 rounded bg-sky-500/10 text-sky-400 border border-sky-500/30">curée</span>`
    : `<span class="text-[10px] px-1.5 py-0.5 rounded bg-slate-700/50 text-slate-400 border border-slate-600">générée</span>`;
}

function newBadge(firstSeen) {
  const today = new Date().toISOString().slice(0, 10);
  if (firstSeen !== today) return "";
  return `<span class="text-[10px] px-1.5 py-0.5 rounded bg-fuchsia-500/10 text-fuchsia-400 border border-fuchsia-500/30">nouveau</span>`;
}

function toxicityBadge(alerts) {
  if (!alerts || alerts.length === 0) return "";
  return `<span title="${alerts.join(', ')}" class="text-[10px] px-1.5 py-0.5 rounded bg-orange-500/10 text-orange-400 border border-orange-500/30">⚠ BRENK ×${alerts.length}</span>`;
}

function moleculeCard(m) {
  return `<article data-id="${m.id}" class="mol-card cursor-pointer bg-slate-900 border border-slate-800 hover:border-emerald-500/50 rounded-xl p-4 transition">
    <div class="flex items-start justify-between gap-2">
      <h3 class="font-mono text-sm text-emerald-300 truncate">${m.id}</h3>
      <div class="flex items-center gap-1 shrink-0">${newBadge(m.first_seen)}${sourceBadge(m.source)}${badge(m.tpp_pass)}</div>
    </div>
    <div class="mol-thumb bg-white rounded-lg h-32 my-3 flex items-center justify-center overflow-hidden" data-smiles="${m.canonical_smiles}"></div>
    <dl class="grid grid-cols-3 gap-x-2 gap-y-1 text-xs text-slate-400">
      <div>MW <span class="text-slate-200 font-mono">${m.mw}</span></div>
      <div>LogP <span class="text-slate-200 font-mono">${m.logp}</span></div>
      <div>TPSA <span class="text-slate-200 font-mono">${m.tpsa}</span></div>
      <div>HBD <span class="text-slate-200 font-mono">${m.hbd}</span></div>
      <div>HBA <span class="text-slate-200 font-mono">${m.hba}</span></div>
      <div>SA <span class="text-slate-200 font-mono">${m.sa_score}</span></div>
    </dl>
    <div class="mt-2">${toxicityBadge(m.toxicity_alerts)}</div>
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

async function drawThumbnail2D(el) {
  try {
    const RDKit = await getRDKit();
    const mol = RDKit.get_mol(el.dataset.smiles);
    el.innerHTML = mol.get_svg(200, 120);
    mol.delete();
  } catch (e) {
    el.innerHTML = `<span class="text-slate-400 text-xs">rendu 2D indisponible</span>`;
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
    el.innerHTML = `<span class="text-slate-400 text-xs">rendu 2D indisponible</span>`;
  }
}

// --- 3Dmol.js : rendu du conformère 3D depuis le bloc SDF précalculé ---
function draw3D(sdf) {
  const container = document.getElementById("modal-3d");
  container.innerHTML = "";
  if (!sdf) {
    container.innerHTML = `<div class="h-full flex items-center justify-center text-slate-500 text-xs">Pas de conformère 3D</div>`;
    return;
  }
  const viewer = $3Dmol.createViewer(container, { backgroundColor: "#020617" });
  viewer.addModel(sdf, "sdf");
  viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });
  viewer.zoomTo();
  viewer.render();
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

async function openModal(id) {
  const m = state.molecules.find((mol) => mol.id === id);
  if (!m) return;
  document.getElementById("modal-title").textContent = m.id;
  document.getElementById("modal-props").innerHTML = [
    ["MW", m.mw], ["LogP", m.logp], ["TPSA", m.tpsa], ["HBD", m.hbd],
    ["HBA", m.hba], ["RotB", m.rotatable_bonds], ["SA score", m.sa_score], ["QED", m.qed],
    ["Fitness", m.fitness ?? "—"], ["Découverte le", m.first_seen ?? "—"],
  ].map(([k, v]) => `<div class="bg-slate-800/60 rounded-lg px-2 py-1"><span class="text-slate-400">${k}:</span> <span class="font-mono">${v}</span></div>`).join("");
  document.getElementById("modal-notes").textContent =
    `${m.formula} — ${m.notes}` +
    (m.toxicity_alerts && m.toxicity_alerts.length ? ` | Alertes BRENK : ${m.toxicity_alerts.join(", ")}` : "");
  document.getElementById("modal").classList.remove("hidden");
  document.getElementById("modal").classList.add("flex");

  drawModal2D(m.canonical_smiles);
  document.getElementById("modal-3d").innerHTML = `<div class="h-full flex items-center justify-center text-slate-500 text-xs">Chargement du conformère…</div>`;
  const conformers = await getConformers();
  draw3D(conformers[m.id]);
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
    const matchesQuery = m.id.toLowerCase().includes(q) || m.smiles.toLowerCase().includes(q);
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
  renderSummary(payload.molecules);
  renderChart(payload.molecules);
  applyFilters();
}

main().catch((err) => {
  console.error(err);
  document.getElementById("molecule-grid").innerHTML =
    `<p class="text-red-400 col-span-full">Erreur de chargement de data/molecules.json — ${err.message}</p>`;
});
