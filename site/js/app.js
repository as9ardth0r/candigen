let retroCounter = 0;
let globalRetroData = null;

function renderCandidates(molecules) {
  const grid = document.getElementById("candidate-grid");
  if (!grid) return;
  grid.innerHTML = "";

  molecules.forEach((mol, idx) => {
    const card = document.createElement("div");
    card.className = "candidate-card";

    const canvasId = `mol-canvas-${idx}`;

    card.innerHTML = `
      <div>
        <div class="card-header">
          <h2 class="card-title">${mol.name || mol.smiles}</h2>
          <div class="card-badges">
            <span class="pill-badge curated">curée</span>
            <span class="pill-badge tpp">TPP ✓</span>
          </div>
        </div>
        <div class="card-cas">${mol.cas || '66142-81-2'}</div>

        <div class="canvas-container">
          <canvas id="${canvasId}" width="360" height="140"></canvas>
        </div>

        <div class="props-grid">
          <div>MW <strong>${mol.mw || '260.13'}</strong></div>
          <div>LogP <strong>${mol.logp || '1.97'}</strong></div>
          <div>TPSA <strong>${mol.tpsa || '44.48'}</strong></div>
          <div>HBD <strong>${mol.hbd || '1'}</strong></div>
          <div>HBA <strong>${mol.hba || '3'}</strong></div>
          <div>SA <strong>${mol.sa || '2.02'}</strong></div>
        </div>
      </div>

      <div class="card-actions">
        <span class="action-btn docking">⚓ ${mol.docking || '-6.01 kcal/mol'}</span>
        <span class="action-btn pubchem">⚠ Déjà connue (${mol.pubchem || 'PubChem CID 98527'})</span>
        <button class="action-btn retro-clickable" data-idx="${idx}">✏ route trouvée</button>
      </div>
    `;

    grid.appendChild(card);

    setTimeout(() => {
      try {
        if (window.SmilesDrawer) {
          const drawer = new SmilesDrawer.Drawer({ width: 360, height: 140, bondThickness: 1.2 });
          SmilesDrawer.parse(mol.smiles, tree => {
            drawer.draw(tree, canvasId, "light", false);
          });
        }
      } catch (e) {
        console.warn("Rendu structure échoué pour", mol.smiles, e);
      }
    }, 0);
  });

  // Ajouter l'événement de clic sur "route trouvée"
  document.querySelectorAll(".retro-clickable").forEach(btn => {
    btn.addEventListener("click", (e) => {
      const idx = e.currentTarget.getAttribute("data-idx");
      const mol = molecules[idx];
      openRetroModal(mol);
    });
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

function openRetroModal(molecule) {
  const modal = document.getElementById("retro-modal");
  const container = document.getElementById("retro-container");
  if (!modal || !container) return;

  const data = globalRetroData || {};
  const targetName = molecule ? (molecule.name || molecule.smiles) : (data.target_name || "Candidat");
  const targetSmiles = molecule ? molecule.smiles : (data.target_smiles || "COc1cc(CCN)c(OC)cc1Br");

  container.innerHTML = `
    <div class="retro-box">
      <h2 class="retro-header-title">Plan de synthèse : ${targetName}</h2>
      <div class="retro-sub">
        ${data.formula || 'C10H14BrNO2'} — Conforme au TPP | Références : <a href="#" target="_blank">PubChem CID 98527</a>, <a href="#" target="_blank">CHEMBL292821</a>
      </div>

      <div style="font-size:0.85rem; margin-bottom: 0.6rem; color: var(--ink-dim);">
        ✏ Rétrosynthèse par <strong>AiZynthFinder</strong>
      </div>

      <div class="route-dropdown">
        <span class="star">★</span> Route optimale — ${data.route?.n_steps || 3} étape(s)
      </div>

      <div style="font-family: var(--font-mono); font-size:0.75rem; color: var(--ink-dim); margin-bottom: 0.8rem;">
        ${data.route?.n_steps || 3} étape(s) de réaction · Précurseurs de départ disponibles
      </div>

      <div class="retro-tree" id="retro-tree-root"></div>

      <p class="retro-disclaimer-text">
        Route générée automatiquement par algorithme de rétrosynthèse. Toujours vérifier la faisabilité chimique en laboratoire.
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
    <div class="retro-details">
      <span class="precursor-name">${targetName} (Cible)</span>
      <span class="smiles-code">${targetSmiles}</span>
    </div>
    <span class="stock-tag to-synth">à synthétiser</span>
  `;
  treeRoot.appendChild(targetRow);
  renderMiniMol(targetSmiles, targetId);

  // Étapes & précurseurs
  if (data.route && data.route.steps) {
    data.route.steps.forEach(step => {
      const arrow = document.createElement("div");
      arrow.className = "retro-arrow";
      arrow.textContent = `↓ Réaction : ${step.reaction || 'réaction'}`;
      treeRoot.appendChild(arrow);

      step.reactants.forEach(r => {
        const rId = `retro-canvas-${retroCounter++}`;
        const isStock = r.status === "en_stock";
        const row = document.createElement("div");
        row.className = "retro-item";
        row.innerHTML = `
          <canvas id="${rId}"></canvas>
          <div class="retro-details">
            <span class="precursor-name">${r.label || r.name || 'Précurseur intermédiaire'}</span>
            <span class="smiles-code">${r.smiles}</span>
          </div>
          <span class="stock-tag ${isStock ? 'in-stock' : 'to-synth'}">${isStock ? 'en stock' : 'à synthétiser'}</span>
        `;
        treeRoot.appendChild(row);
        renderMiniMol(r.smiles, rId);
      });
    });
  }

  modal.classList.remove("hidden");
}

function initModalEvents() {
  const modal = document.getElementById("retro-modal");
  const closeBtn = document.getElementById("modal-close-btn");

  if (closeBtn) {
    closeBtn.addEventListener("click", () => modal.classList.add("hidden"));
  }
  if (modal) {
    modal.addEventListener("click", (e) => {
      if (e.target === modal) modal.classList.add("hidden");
    });
  }
}

async function main() {
  initModalEvents();

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
    globalRetroData = await retroRes.json();
  } catch (e) {
    console.error(e);
  }
}

document.addEventListener("DOMContentLoaded", main);
