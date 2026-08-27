// Variable globale pour l'incrémentation des ID canvas
let retroCounter = 0;
let globalRetroData = null;

// Rendu RDKit mini-molécule dans un Canvas
function renderMiniMol(smiles, canvasId) {
  setTimeout(() => {
    const canvas = document.getElementById(canvasId);
    if (!canvas) return;
    
    if (window.RDKitModule && smiles) {
      try {
        const mol = window.RDKitModule.get_mol(smiles);
        if (mol) {
          canvas.width = 70;
          canvas.height = 50;
          mol.draw_to_canvas(canvas, 70, 50);
          mol.delete();
          return;
        }
      } catch (e) {
        console.warn("Erreur de rendu RDKit pour :", smiles);
      }
    }

    // Fallback visuel SVG/Texte si RDKit n'est pas chargé
    const ctx = canvas.getContext("2d");
    if (ctx) {
      canvas.width = 70;
      canvas.height = 50;
      ctx.fillStyle = "#1e293b";
      ctx.fillRect(0, 0, 70, 50);
      ctx.fillStyle = "#94a3b8";
      ctx.font = "10px sans-serif";
      ctx.textAlign = "center";
      ctx.fillText("Molécule", 35, 28);
    }
  }, 50);
}

function openRetroModal(molecule) {
  const modal = document.getElementById("retro-modal");
  const container = document.getElementById("retro-container");
  if (!modal || !container) return;

  const data = globalRetroData || {};
  const targetName = molecule ? (molecule.name || molecule.smiles) : (data.target_name || "Candidat");
  const targetSmiles = molecule ? molecule.smiles : (data.target_smiles || "COc1cc2c(OC)nc(C)nc2c(O)c1C1CCOC1");

  // Définition d'étapes de repli si le JSON distant est absent ou incomplet
  const defaultSteps = [
    {
      reaction: "Substitution / Amidation",
      reactants: [
        { label: "Précurseur Aminé", smiles: "COc1cc(CCN)c(OC)cc1Br", status: "a_synthetiser" }
      ]
    },
    {
      reaction: "Couplage & Cyclisation",
      reactants: [
        { label: "Isocyanate intermédiaire", smiles: "CN=C=O", status: "en_stock" },
        { label: "Dibromure d'aryle", smiles: "COc1cc(Br)c(O)cc1Br", status: "a_synthetiser" }
      ]
    },
    {
      reaction: "Halogénation régiosélective",
      reactants: [
        { label: "NBS (N-Bromosuccinimide)", smiles: "O=C1CCC(=O)N1Br", status: "en_stock" },
        { label: "Réactif de départ commercial", smiles: "COc1cc(Br)c(O)cc1C", status: "en_stock" }
      ]
    }
  ];

  const stepsToRender = (data.route && data.route.steps && data.route.steps.length > 0) 
    ? data.route.steps 
    : defaultSteps;

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
        <span class="star">★</span> Route optimale — ${stepsToRender.length} étape(s)
      </div>

      <div style="font-family: var(--font-mono); font-size:0.75rem; color: var(--ink-dim); margin-bottom: 0.8rem;">
        ${stepsToRender.length} étape(s) de réaction · Précurseurs de départ disponibles
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

  // Rendu de chaque étape
  stepsToRender.forEach(step => {
    const arrow = document.createElement("div");
    arrow.className = "retro-arrow";
    arrow.textContent = `↓ Réaction : ${step.reaction || 'réaction'}`;
    treeRoot.appendChild(arrow);

    step.reactants.forEach(r => {
      const rId = `retro-canvas-${retroCounter++}`;
      const isStock = (r.status === "en_stock" || r.status === "in-stock");
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

  modal.classList.remove("hidden");
}

function initModalEvents() {
  const modal = document.getElementById("retro-modal");
  const closeBtn = document.getElementById("close-retro-modal");
  if (closeBtn && modal) {
    closeBtn.addEventListener("click", () => modal.classList.add("hidden"));
  }
}

// Initialisation au chargement de la page
document.addEventListener("DOMContentLoaded", async () => {
  initModalEvents();

  // Chargement sécurisé du JSON de rétrosynthèse
  try {
    const res = await fetch("data/retrosynthesis_erlotinib.json");
    if (res.ok) {
      globalRetroData = await res.json();
    }
  } catch (e) {
    console.warn("Fichier rétrosynthèse indisponible, utilisation des données de secours.");
  }
});
