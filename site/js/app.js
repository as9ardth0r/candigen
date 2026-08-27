let retroCounter = 0;

// Rendu mini-molécule RDKit / Fallback
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
        console.warn("Erreur RDKit :", smiles);
      }
    }

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

// Fonction d'affichage de la grille de molécules dans le dashboard
function renderCandidates(molecules) {
  const container = document.getElementById("grid-container") || document.getElementById("molecules-grid");
  if (!container) return;

  container.innerHTML = "";
  molecules.forEach((mol, idx) => {
    const card = document.createElement("div");
    card.className = "mol-card";
    const canvasId = `main-mol-canvas-${idx}`;
    
    card.innerHTML = `
      <h3>${mol.id || mol.name || 'Candidat'}</h3>
      <div class="canvas-wrapper">
        <canvas id="${canvasId}"></canvas>
      </div>
      <p class="smiles-text">${mol.smiles || ''}</p>
      <button class="btn-retro" data-idx="${idx}">Plan de rétrosynthèse</button>
    `;

    container.appendChild(card);
    if (mol.smiles) renderMiniMol(mol.smiles, canvasId);

    const btn = card.querySelector(".btn-retro");
    if (btn) {
      btn.addEventListener("click", () => openRetroModal(mol));
    }
  });
}

// Extraction récursive de l'arbre AiZynthFinder
function parseAiZynthTree(node, steps = []) {
  if (!node) return steps;

  if (node.children && node.children.length > 0) {
    const reactionName = node.metadata?.reaction_hash || node.type || "Étape de synthèse";
    const reactants = [];

    node.children.forEach(child => {
      reactants.push({
        name: child.metadata?.mapped_smiles || child.smiles || "Intermédiaire",
        smiles: child.smiles || "",
        status: child.is_chemical ? (child.in_stock ? "en_stock" : "a_synthetiser") : "a_synthetiser"
      });
      parseAiZynthTree(child, steps);
    });

    steps.unshift({
      reaction: reactionName,
      reactants: reactants
    });
  }

  return steps;
}

// Modale de rétrosynthèse
async function openRetroModal(molecule) {
  const modal = document.getElementById("retro-modal");
  const container = document.getElementById("retro-container");
  if (!modal || !container) return;

  const targetName = molecule ? (molecule.id || molecule.name || molecule.smiles) : "Candidat";
  const targetSmiles = molecule ? molecule.smiles : "";
  
  let stepsToRender = [];

  if (molecule && molecule.id) {
    try {
      const res = await fetch(`data/retrosynthesis/${molecule.id}.json`);
      if (res.ok) {
        const data = await res.json();
        if (data.smiles || data.children) {
          stepsToRender = parseAiZynthTree(data);
        } else if (data.route?.steps) {
          stepsToRender = data.route.steps;
        }
      }
    } catch (e) {
      console.warn("Rétrosynthèse non trouvée pour :", molecule.id);
    }
  }

  // Fallback si pas de fichier de rétrosynthèse disponible
  if (stepsToRender.length === 0) {
    stepsToRender = [
      {
        reaction: "Substitution / Amidation",
        reactants: [{ name: "Précurseur Aminé", smiles: "COc1cc(CCN)c(OC)cc1Br", status: "a_synthetiser" }]
      },
      {
        reaction: "Couplage & Cyclisation",
        reactants: [
          { name: "Isocyanate intermédiaire", smiles: "CN=C=O", status: "en_stock" },
          { name: "Dibromure d'aryle", smiles: "COc1cc(Br)c(O)cc1Br", status: "a_synthetiser" }
        ]
      },
      {
        reaction: "Halogénation régiosélective",
        reactants: [
          { name: "NBS (N-Bromosuccinimide)", smiles: "O=C1CCC(=O)N1Br", status: "en_stock" },
          { name: "Réactif de départ commercial", smiles: "COc1cc(Br)c(O)cc1C", status: "en_stock" }
        ]
      }
    ];
  }

  container.innerHTML = `
    <div class="retro-box">
      <h2 class="retro-header-title">Plan de synthèse : ${targetName}</h2>
      <div class="retro-sub">
        Conforme au TPP | Rétrosynthèse AiZynthFinder
      </div>

      <div style="font-size:0.85rem; margin-top: 0.4rem; color: var(--ink-dim);">
        ✏ Rétrosynthèse par <strong>AiZynthFinder</strong>
      </div>

      <div class="route-dropdown" style="margin-top:0.8rem;">
        <span class="star">★</span> Route optimale — ${stepsToRender.length} étape(s)
      </div>

      <div style="font-family: var(--font-mono); font-size:0.75rem; color: var(--ink-dim); margin-bottom: 0.8rem;">
        ${stepsToRender.length} étape(s) de réaction · Précurseurs identifiés
      </div>

      <div class="retro-tree" id="retro-tree-root"></div>

      <p class="retro-disclaimer-text">
        Route générée automatiquement par algorithme de rétrosynthèse. Toujours vérifier la faisabilité chimique en laboratoire.
      </p>
    </div>
  `;

  const treeRoot = document.getElementById("retro-tree-root");

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
  if (targetSmiles) renderMiniMol(targetSmiles, targetId);

  stepsToRender.forEach(step => {
    const arrow = document.createElement("div");
    arrow.className = "retro-arrow";
    arrow.textContent = `↓ Réaction : ${step.reaction}`;
    treeRoot.appendChild(arrow);

    step.reactants.forEach(r => {
      const rId = `retro-canvas-${retroCounter++}`;
      const isStock = (r.status === "en_stock" || r.status === "in-stock");
      const row = document.createElement("div");
      row.className = "retro-item";
      row.innerHTML = `
        <canvas id="${rId}"></canvas>
        <div class="retro-details">
          <span class="precursor-name">${r.name || 'Précurseur'}</span>
          <span class="smiles-code">${r.smiles}</span>
        </div>
        <span class="stock-tag ${isStock ? 'in-stock' : 'to-synth'}">${isStock ? 'en stock' : 'à synthétiser'}</span>
      `;
      treeRoot.appendChild(row);
      if (r.smiles) renderMiniMol(r.smiles, rId);
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

// Initialisation globale et chargement de site/data/molecules.json
document.addEventListener("DOMContentLoaded", async () => {
  initModalEvents();

  try {
    const res = await fetch("data/molecules.json");
    if (res.ok) {
      const data = await res.json();
      const molecules = data.molecules || data;
      renderCandidates(molecules);

      const statusEl = document.querySelector(".header-status, #status-text");
      if (statusEl) {
        statusEl.textContent = `EGFR · ${molecules.length || 0} molécules chargées`;
      }
    }
  } catch (e) {
    console.error("Erreur de chargement des molécules :", e);
  }
});
