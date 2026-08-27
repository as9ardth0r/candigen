const HEAT_LOW = [59, 130, 196];   // bleu
const HEAT_MID = [232, 180, 76];   // ambre
const HEAT_HIGH = [228, 87, 46];   // rouge-orangé

const prefersReducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

function heatColor(t) {
  // t dans [0, 1] -> interpolation bleu -> ambre -> rouge (façon facteur B)
  t = Math.max(0, Math.min(1, t));
  const [a, b] = t < 0.5 ? [HEAT_LOW, HEAT_MID] : [HEAT_MID, HEAT_HIGH];
  const local = t < 0.5 ? t / 0.5 : (t - 0.5) / 0.5;
  const rgb = a.map((c, i) => Math.round(c + (b[i] - c) * local));
  return `rgb(${rgb.join(",")})`;
}

function renderStats(data) {
  if (!data || !data.molecules) return;
  const best = Math.max(...data.molecules.map(m => m.fitness), 0);
  
  const elGen = document.getElementById("stat-generations");
  const elBest = document.getElementById("stat-best");
  const elPop = document.getElementById("stat-population");

  if (elGen) elGen.textContent = data.n_generations;
  if (elBest) elBest.textContent = best.toFixed(3);
  if (elPop) elPop.textContent = data.molecules.length;
}

function renderChart(history) {
  const container = document.getElementById("chart-container");
  if (!container) return;

  if (!history || history.length === 0) {
    container.innerHTML = '<p style="color:#7C8798">Pas d\'historique disponible.</p>';
    return;
  }

  const width = 880, height = 260, padding = 36;
  const maxGen = history.length - 1;

  const x = i => padding + (i / Math.max(1, maxGen)) * (width - 2 * padding);
  const y = v => height - padding - v * (height - 2 * padding);

  const bestLine = history.map((h, i) => `${x(i)},${y(h.best_fitness)}`).join(" L ");
  const meanPoints = history.map((h, i) => `${x(i)},${y(h.mean_fitness)}`).join(" ");
  const lastBest = history[history.length - 1].best_fitness;
  const lastColor = heatColor(lastBest);

  // aire sous la courbe "meilleure fitness", refermée sur la ligne de base
  const areaPath = `M ${x(0)},${y(0)} L ${bestLine} L ${x(maxGen)},${y(0)} Z`;
  const lineDrawClass = prefersReducedMotion ? "" : "chart-line-draw";

  const gridLines = [0, 0.25, 0.5, 0.75, 1].map(v => `
    <line x1="${padding}" y1="${y(v)}" x2="${width - padding}" y2="${y(v)}"
          stroke="rgba(148,163,184,0.14)" stroke-width="1" />
    <text x="${padding - 8}" y="${y(v) + 4}" text-anchor="end"
          font-family="IBM Plex Mono" font-size="10" fill="#7C8798">${v.toFixed(2)}</text>
  `).join("");

  container.innerHTML = `
    <svg viewBox="0 0 ${width} ${height}">
      <defs>
        <linearGradient id="areaGradient" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0%" stop-color="${lastColor}" stop-opacity="0.28" />
          <stop offset="100%" stop-color="${lastColor}" stop-opacity="0" />
        </linearGradient>
        <linearGradient id="lineGradient" x1="0" y1="0" x2="1" y2="0">
          <stop offset="0%" stop-color="rgb(${HEAT_LOW.join(",")})" />
          <stop offset="50%" stop-color="rgb(${HEAT_MID.join(",")})" />
          <stop offset="100%" stop-color="rgb(${HEAT_HIGH.join(",")})" />
        </linearGradient>
        <filter id="glow" x="-40%" y="-40%" width="180%" height="180%">
          <feGaussianBlur stdDeviation="2.4" result="blur" />
          <feMerge>
            <feMergeNode in="blur" />
            <feMergeNode in="SourceGraphic" />
          </feMerge>
        </filter>
      </defs>

      ${gridLines}
      <path d="${areaPath}" fill="url(#areaGradient)" stroke="none" />
      <polyline points="${meanPoints}" fill="none" stroke="#3A4152" stroke-width="2" stroke-dasharray="4 4" />
      <path d="M ${bestLine}" fill="none" stroke="url(#lineGradient)" stroke-width="2.5"
            filter="url(#glow)" pathLength="1" class="${lineDrawClass}" />
      ${history.map((h, i) => `<circle cx="${x(i)}" cy="${y(h.best_fitness)}" r="3" fill="${heatColor(h.best_fitness)}" />`).join("")}
    </svg>
    <div class="chart-legend">
      <span><span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:${lastColor};box-shadow:0 0 6px ${lastColor};margin-right:5px"></span>meilleure fitness</span>
      <span><span style="display:inline-block;width:10px;height:2px;background:#3A4152;margin-right:5px;vertical-align:middle"></span>fitness moyenne</span>
    </div>
  `;
}

function renderCandidates(molecules) {
  const grid = document.getElementById("candidate-grid");
  if (!grid) return;
  grid.innerHTML = "";

  molecules.forEach((mol, idx) => {
    const card = document.createElement("div");
    card.className = "candidate-card";

    const rankTag = document.createElement("span");
    rankTag.className = "rank-tag";
    rankTag.textContent = `#${String(idx + 1).padStart(2, "0")}`;

    const canvas = document.createElement("canvas");
    canvas.id = `mol-canvas-${idx}`;
    canvas.width = 300;
    canvas.height = 200;

    const smilesEl = document.createElement("p");
    smilesEl.className = "candidate-smiles";
    smilesEl.textContent = mol.smiles;

    const badgeRow = document.createElement("div");
    badgeRow.className = "badge-row";
    const fitnessBadge = document.createElement("span");
    fitnessBadge.className = "fitness-badge";
    fitnessBadge.textContent = mol.fitness.toFixed(3);
    const color = heatColor(mol.fitness);
    fitnessBadge.style.background = color;
    fitnessBadge.style.boxShadow = `0 0 10px 0 ${color.replace("rgb", "rgba").replace(")", ",0.55)")}`;
    const qedSpan = document.createElement("span");
    qedSpan.style.fontFamily = "var(--font-mono)";
    qedSpan.style.fontSize = "0.78rem";
    qedSpan.style.color = "var(--ink-dim)";
    qedSpan.textContent = `QED ${mol.qed.toFixed(2)}`;
    badgeRow.append(fitnessBadge, qedSpan);

    const tagRow = document.createElement("div");
    tagRow.className = "tag-row";
    tagRow.innerHTML = `
      <span class="tag ${mol.drug_like ? "ok" : "warn"}">${mol.drug_like ? "drug-like" : "hors Lipinski"}</span>
      <span class="tag ${mol.structurally_clean ? "ok" : "warn"}">${mol.structurally_clean ? "propre" : "alerte structurale"}</span>
    `;

    card.append(rankTag, canvas, smilesEl, badgeRow, tagRow);
    grid.appendChild(card);

    try {
      if (window.SmilesDrawer) {
        const drawer = new SmilesDrawer.Drawer({ width: 300, height: 200, bondThickness: 1.2 });
        SmilesDrawer.parse(mol.smiles, tree => {
          drawer.draw(tree, canvas.id, "light", false);
        }, err => console.warn("SmilesDrawer parse error", mol.smiles, err));
      }
    } catch (e) {
      console.warn("Rendu structure échoué pour", mol.smiles, e);
    }
  });
}

let retroMolCounter = 0;

function renderMolCard(mol, container) {
  if (!container) return;
  const card = document.createElement("div");
  card.className = "mol-card";

  const canvasId = `retro-mol-${retroMolCounter++}`;
  const canvas = document.createElement("canvas");
  canvas.id = canvasId;
  canvas.width = 220;
  canvas.height = 150;

  const label = document.createElement("p");
  label.className = "mol-label";
  label.textContent = mol.label;

  card.append(canvas, label);

  if (mol.status) {
    const badge = document.createElement("span");
    const isStock = mol.status === "en_stock";
    badge.className = `status-badge ${isStock ? "en-stock" : "a-synthetiser"}`;
    badge.textContent = isStock ? "en stock" : "à synthétiser";
    card.appendChild(badge);
  }

  container.appendChild(card);

  try {
    if (window.SmilesDrawer) {
      const drawer = new SmilesDrawer.Drawer({ width: 220, height: 150, bondThickness: 1.1 });
      SmilesDrawer.parse(mol.smiles, tree => {
        drawer.draw(tree, canvasId, "light", false);
      }, err => console.warn("SmilesDrawer parse error", mol.smiles, err));
    }
  } catch (e) {
    console.warn("Rendu structure échoué pour", mol.smiles, e);
  }
}

function renderRetrosynthesis(data) {
  const container = document.getElementById("retro-container");
  if (!container) return;
  container.innerHTML = "";

  const header = document.createElement("div");
  header.className = "retro-header";
  header.innerHTML = `
    <div>
      <p class="retro-title">${data.target_name} — ${data.formula} (MW ${data.molecular_weight})</p>
      <p class="retro-meta">${data.known_compound?.note ?? ""}</p>
    </div>
    <div class="route-selector"><span class="star">★</span> ${data.route.label}</div>
  `;
  container.appendChild(header);

  const summary = document.createElement("p");
  summary.className = "route-summary";
  summary.textContent = `${data.route.n_steps} étape(s)`;
  container.appendChild(summary);

  const chain = document.createElement("div");
  chain.className = "route-chain";
  container.appendChild(chain);

  // molécule cible en tête de chaîne
  const targetWrap = document.createElement("div");
  targetWrap.className = "step-molecules";
  renderMolCard({ smiles: data.target_smiles, label: data.target_name, status: null }, targetWrap);
  chain.appendChild(targetWrap);

  data.route.steps.forEach(step => {
    const arrow = document.createElement("div");
    arrow.className = "reaction-arrow";
    arrow.textContent = `réaction — ${step.reaction}`;
    chain.appendChild(arrow);

    const stepWrap = document.createElement("div");
    stepWrap.className = "step-molecules";
    step.reactants.forEach(r => renderMolCard(r, stepWrap));
    chain.appendChild(stepWrap);
  });

  const disclaimer = document.createElement("p");
  disclaimer.className = "retro-disclaimer";
  disclaimer.textContent = data.disclaimer;
  container.appendChild(disclaimer);
}

async function main() {
  try {
    const res = await fetch("data/molecules.json");
    const data = await res.json();

    renderStats(data);
    renderChart(data.history);
    renderCandidates(data.molecules);

    const date = new Date(data.generated_at);
    const metaEl = document.getElementById("footer-meta");
    if (metaEl) {
      metaEl.textContent = `Cible : ${data.target} — généré le ${date.toLocaleDateString("fr-FR")} à ${date.toLocaleTimeString("fr-FR")}`;
    }
  } catch (e) {
    const metaEl = document.getElementById("footer-meta");
    if (metaEl) {
      metaEl.textContent = "Aucune donnée disponible — lancez scripts/run_pipeline.py pour générer site/data/molecules.json.";
    }
    console.error(e);
  }

  try {
    const retroRes = await fetch("data/retrosynthesis_erlotinib.json");
    const retroData = await retroRes.json();
    renderRetrosynthesis(retroData);
  } catch (e) {
    const container = document.getElementById("retro-container");
    if (container) {
      container.innerHTML = '<p style="color:var(--ink-dim)">Pas de route disponible.</p>';
    }
    console.error(e);
  }
}

document.addEventListener("DOMContentLoaded", main);
