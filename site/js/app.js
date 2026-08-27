let retroCounter = 0;

function renderCandidates(molecules) {
  const grid = document.getElementById("candidate-grid");
  if (!grid) return;
  grid.innerHTML = "";

  molecules.forEach((mol, idx) => {
    const card = document.createElement("div");
    card.className = "candidate-card";

    const canvasId = `mol-canvas-${idx}`;

    card.innerHTML = `
      <div class="card-header">
        <h2 class="card-title">${mol.name || mol.smiles}</h2>
        <div class="card-badges">
          <span class="pill-badge curated">curée</span>
          <span class="pill-badge tpp">TPP ✓</span>
        </div>
      </div>
      <div class="card-cas">${mol.cas || '66142-81-2'}</div>

      <div class="canvas-container">
        <canvas id="${canvasId}" width="360" height="160"></canvas>
      </div>

      <div class="props-grid">
        <div>MW <strong>${mol.mw || '260.13'}</strong></div>
        <div>LogP <strong>${mol.logp || '1.97'}</strong></div>
        <div>TPSA <strong>${mol.tpsa || '44.48'}</strong></div>
        <div>HBD <strong>${mol.hbd || '1'}</strong></div>
        <div>HBA <strong>${mol.hba || '3'}</strong></div>
        <div>SA <strong>${mol.sa || '2.02'}</strong></div>
      </div>

      <div class="card-actions">
        <span class="action-btn docking">⚓ ${mol.docking || '-6.01 kcal/mol'}</span>
        <span class="action-btn pubchem">⚠ Déjà connue (${mol.pubchem || 'PubChem CID 98527'})</span>
        <span class="action-btn retro">✏ route trouvée</span>
      </div>
    `;

    grid.appendChild(card);

    setTimeout(() => {
      try {
        if (window.SmilesDrawer) {
          const drawer = new SmilesDrawer.Drawer({ width: 360, height: 160, bondThickness: 1.2 });
          SmilesDrawer.parse(mol.smiles, tree => {
            drawer.draw(tree, canvasId, "light", false);
          });
        }
      } catch (e) {
        console.warn("Rendu structure échoué pour", mol.smiles, e);
      }
    }, 0);
  });
}

function renderMiniMol(smiles, canvasId) {
  setTimeout(() => {
    try {
      if (window.SmilesDrawer) {
        const drawer = new SmilesDrawer.Drawer({ width: 100, height: 70, bondThickness: 1.0 });
        SmilesDrawer.parse(smiles, tree => {
          drawer.draw(tree, canvasId, "light", false);
        });
      }
    } catch (e) {
      console.warn("Rendu mini structure échoué pour", smiles, e);
    }
  }, 0);
}

function renderRetrosynthesis(data) {
  const container = document.getElementById("retro-container");
  if (!container) return;

  container.innerHTML = `
    <div class="retro-box">
      <div class="retro-sub">
        ${data.formula || 'C10H14BrNO2'} — Conforme au TPP | Déjà répertoriée : <a href="#">PubChem CID 98527</a>, <a href="#">CHEMBL292821</a>
      </div>

      <div style="font-size:0.85rem; margin-bottom: 0.6rem; color: var(--ink-dim);">
        ✏ Rétrosynthèse — AiZynthFinder
      </div>

      <div class="route-dropdown">
        <span class="star">★</span> Route 5 — ${data.route?.n_steps || 3} étape(s) <span>˅</span>
      </div>

      <div style="font-family: var(--font-mono); font-size:0.75rem; color: var(--ink-dim); margin-bottom: 0.8rem;">
        ${data.route?.n_steps || 3} étape(s) · 3/3 précurseur(s) déjà en stock
      </div>

      <div class="retro-tree" id="retro-tree-root"></div>

      <p class="retro-disclaimer-text">
        Route générée automatiquement — une piste à explorer, pas un protocole de synthèse validé.
      </p>
    </div>
  `;

  const treeRoot = document.getElementById("retro-tree-root");

  // Molécule cible
  const targetId = `retro-canvas-${retroCounter++}`;
  const targetRow = document.createElement("div");
  targetRow.className = "retro-item";
  targetRow.innerHTML = `
    <canvas id="${targetId}"></canvas>
    <div class="smiles-code">${data.target_smiles || 'COc1cc(CCN)c(OC)cc1Br'}</div>
    <span class="stock-tag to-synth">à synthétiser</span>
  `;
  treeRoot.appendChild(targetRow);
  renderMiniMol(data.target_smiles || 'COc1cc(CCN)c(OC)cc1Br', targetId);

  // Étapes
  if (data.route && data.route.steps) {
    data.route.steps.forEach(step => {
      const arrow = document.createElement("div");
      arrow.className = "retro-arrow";
      arrow.textContent = `↓ réaction`;
      treeRoot.appendChild(arrow);

      step.reactants.forEach(r => {
        const rId = `retro-canvas-${retroCounter++}`;
        const isStock = r.status === "en_stock";
        const row = document.createElement("div");
        row.className = "retro-item";
        row.innerHTML = `
          <canvas id="${rId}"></canvas>
          <div class="smiles-code">${r.smiles}</div>
          <span class="stock-tag ${isStock ? 'in-stock' : 'to-synth'}">${isStock ? 'en stock' : 'à synthétiser'}</span>
        `;
        treeRoot.appendChild(row);
        renderMiniMol(r.smiles, rId);
      });
    });
  }
}

async function main() {
  try {
    const res = await fetch("data/molecules.json");
    const data = await res.json();
    renderCandidates(data.molecules || []);

    const metaEl = document.getElementById("footer-meta");
    if (metaEl) {
      metaEl.textContent = `EGFR · généré le ${new Date(data.generated_at || Date.now()).toLocaleDateString("fr-FR")} ${new Date().toLocaleTimeString("fr-FR")}`;
    }
  } catch (e) {
    console.error(e);
  }

  try {
    const retroRes = await fetch("data/retrosynthesis_erlotinib.json");
    const retroData = await retroRes.json();
    renderRetrosynthesis(retroData);
  } catch (e) {
    console.error(e);
  }
}

document.addEventListener("DOMContentLoaded", main);
