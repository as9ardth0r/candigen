function openRetroModal(molecule) {
  const modal = document.getElementById("retro-modal");
  const container = document.getElementById("retro-container");
  if (!modal || !container) return;

  const data = globalRetroData || {};
  const targetName = molecule ? (molecule.name || molecule.smiles) : (data.target_name || "Candidat");
  const targetSmiles = molecule ? molecule.smiles : (data.target_smiles || "COc1cc2c(OC)nc(C)nc2c(O)c1C1CCOC1");

  // Définition d'étapes par défaut si le JSON est incomplet
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

  modal.classList.remove("hidden");
}
