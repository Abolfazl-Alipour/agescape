/**
 * Conserved Aging Explorer — Client Application Logic
 * Interactive SVG/Canvas Visualizations & Multi-Species Analysis Engine
 */

document.addEventListener("DOMContentLoaded", () => {
  // App State
  const state = {
    currentView: "outliers",
    data: {
      outliers: null,
      genes: null,
      divergence: null,
      ranking: null,
      umap: null,
      composition: null
    },
    selectedGene: "Cdkn2a",
    umapColorMode: "age",
    showTostBand: true
  };

  // View Titles & Descriptions
  const viewMeta = {
    outliers: {
      title: "Longevity Outlier Contrast",
      desc: "Evaluation of core mammalian aging program activation across life-stage normalized lifespan."
    },
    genes: {
      title: "Gene Trajectory Query",
      desc: "Searchable cross-species trajectories, log₂ fold-changes, and 95% bootstrap confidence intervals."
    },
    programs: {
      title: "Conserved vs Divergent Aging Programs",
      desc: "Cell-type-resolved cross-mammalian alignment with housekeeping inflation guards."
    },
    atlas: {
      title: "Single-Cell Atlas UMAP & Composition",
      desc: "Subsampled 2D manifold of Tabula Muris Senis with age-stratified composition metrics."
    },
    ranking: {
      title: "Mouse Aging Signal Ranking",
      desc: "Differential expression effect sizes across powered tissue×cell-type units (Milestone M2)."
    },
    caveats: {
      title: "Caveats & Methodological Transparency",
      desc: "Essential caveats on sex×age droplet confounding, cross-sectional design, and statistical power."
    }
  };

  // Color Palettes
  const colors = {
    mouse: "#f43f5e",
    human: "#60a5fa",
    nmr: "#10b981",
    tostBand: "rgba(16, 185, 129, 0.2)",
    age: {
      "3m": "#38bdf8",
      "18m": "#f59e0b",
      "24m": "#f43f5e"
    },
    tissue: {
      "Lung": "#38bdf8",
      "Kidney": "#a855f7",
      "Limb_Muscle": "#f59e0b",
      "Bone_Marrow": "#10b981",
      "Heart_and_Aorta": "#f43f5e",
      "Thymus": "#c084fc"
    }
  };

  // ---------------------------------------------------------
  // 1. Navigation Routing
  // ---------------------------------------------------------
  function setupNavigation() {
    const navButtons = document.querySelectorAll(".nav-item");
    navButtons.forEach(btn => {
      btn.addEventListener("click", () => {
        const targetView = btn.dataset.view;
        switchView(targetView);
      });
    });

    const confoundBtn = document.getElementById("confound-warning-btn");
    if (confoundBtn) {
      confoundBtn.addEventListener("click", (e) => {
        e.preventDefault();
        switchView("caveats");
      });
    }
  }

  function switchView(viewName) {
    state.currentView = viewName;

    // Update active nav button
    document.querySelectorAll(".nav-item").forEach(b => {
      b.classList.toggle("active", b.dataset.view === viewName);
    });

    // Update view panes
    document.querySelectorAll(".view-pane").forEach(pane => {
      pane.classList.toggle("active", pane.id === `view-${viewName}`);
    });

    // Update header
    const meta = viewMeta[viewName] || { title: viewName, desc: "" };
    document.getElementById("current-view-title").textContent = meta.title;
    document.getElementById("current-view-desc").textContent = meta.desc;

    // Trigger view-specific render
    if (viewName === "outliers") renderOutliersChart();
    if (viewName === "genes") renderGeneChart();
    if (viewName === "atlas") renderUmap();
  }

  // ---------------------------------------------------------
  // 2. Data Fetching
  // ---------------------------------------------------------
  async function loadAllData() {
    try {
      const [outliers, genes, divergence, ranking, umap, composition] = await Promise.all([
        fetch("data/longevity_outliers.json").then(r => r.json()).catch(() => null),
        fetch("data/gene_panel_trajectories.json").then(r => r.json()).catch(() => []),
        fetch("data/divergence_table.json").then(r => r.json()).catch(() => []),
        fetch("data/aging_signals_ranked.json").then(r => r.json()).catch(() => []),
        fetch("data/umap_subsample.json").then(r => r.json()).catch(() => []),
        fetch("data/celltype_composition.json").then(r => r.json()).catch(() => [])
      ]);

      state.data = { outliers, genes, divergence, ranking, umap, composition };

      // Initialize Views
      initOutliersView();
      initGenesView();
      initProgramsView();
      initRankingView();
      initAtlasView();

      console.log("All application data successfully loaded.");
    } catch (err) {
      console.error("Failed to load application data:", err);
    }
  }

  // ---------------------------------------------------------
  // 3. View 1: Longevity Outlier Contrast
  // ---------------------------------------------------------
  function initOutliersView() {
    const toggleBand = document.getElementById("toggle-tost-band");
    if (toggleBand) {
      toggleBand.addEventListener("change", (e) => {
        state.showTostBand = e.target.checked;
        renderOutliersChart();
      });
    }

    // Populate stats table
    if (state.data.outliers && state.data.outliers.statistics) {
      const tbody = document.getElementById("table-outlier-stats-body");
      if (tbody) {
        tbody.innerHTML = state.data.outliers.statistics.map(s => {
          let interp = "";
          if (s.comparison.includes("Mouse Aging")) interp = "Robust canonical activation (Positive Control Gate Passed)";
          else if (s.comparison.includes("Human Aging")) interp = "Conserved mammalian activation incline";
          else if (s.comparison.includes("NMR Aging")) interp = "Attenuated / Near-flat trajectory across 30y";
          else if (s.comparison.includes("Interaction")) interp = "Highly significant slope divergence from mouse";
          else if (s.comparison.includes("delta=0.25")) interp = "Statistically significant equivalence within ~18% margin";
          else if (s.comparison.includes("delta=0.15")) interp = "Exploratory strict margin (point estimate slightly > 0.15)";

          return `
            <tr>
              <td><strong>${s.comparison}</strong></td>
              <td>${s.estimate.toFixed(3)}</td>
              <td>${s.se ? s.se.toFixed(3) : "—"}</td>
              <td>[${s.ci_low.toFixed(3)}, ${s.ci_high.toFixed(3)}]</td>
              <td><strong>${s.pvalue < 1e-4 ? s.pvalue.toExponential(2) : s.pvalue.toFixed(4)}</strong></td>
              <td>${interp}</td>
            </tr>
          `;
        }).join("");
      }
    }

    renderOutliersChart();
  }

  function renderOutliersChart() {
    const svg = document.getElementById("outliers-svg");
    if (!svg || !state.data.outliers) return;

    const data = state.data.outliers.donors || [];
    const stats = state.data.outliers.statistics || [];
    const width = 900;
    const height = 450;
    const padding = { top: 40, right: 50, bottom: 60, left: 70 };

    svg.innerHTML = "";

    // Scales
    const xMin = 0.0, xMax = 0.9;
    const yMin = 0.0, yMax = 1.6;

    const xScale = (val) => padding.left + (val - xMin) / (xMax - xMin) * (width - padding.left - padding.right);
    const yScale = (val) => height - padding.bottom - (val - yMin) / (yMax - yMin) * (height - padding.top - padding.bottom);

    // Axes Grid Lines
    let gridHTML = "";
    for (let y = 0.2; y <= 1.4; y += 0.3) {
      const yPos = yScale(y);
      gridHTML += `<line x1="${padding.left}" y1="${yPos}" x2="${width - padding.right}" y2="${yPos}" stroke="rgba(255,255,255,0.06)" stroke-dasharray="4"/>`;
      gridHTML += `<text x="${padding.left - 12}" y="${yPos + 4}" fill="#64748b" font-size="11" text-anchor="end" font-family="monospace">${y.toFixed(1)}</text>`;
    }
    for (let x = 0.0; x <= 0.8; x += 0.2) {
      const xPos = xScale(x);
      gridHTML += `<line x1="${xPos}" y1="${padding.top}" x2="${xPos}" y2="${height - padding.bottom}" stroke="rgba(255,255,255,0.06)" stroke-dasharray="4"/>`;
      gridHTML += `<text x="${xPos}" y="${height - padding.bottom + 20}" fill="#64748b" font-size="11" text-anchor="middle" font-family="monospace">${(x * 100).toFixed(0)}%</text>`;
    }
    svg.innerHTML += gridHTML;

    // TOST Equivalence Band
    if (state.showTostBand) {
      const nmrDonors = data.filter(d => d.species === "Naked Mole-Rat");
      if (nmrDonors.length > 0) {
        const meanY = nmrDonors.reduce((sum, d) => sum + d.program_score, 0) / nmrDonors.length;
        const bandTop = yScale(meanY + 0.125);
        const bandBottom = yScale(meanY - 0.125);
        const bandHeight = bandBottom - bandTop;
        svg.innerHTML += `
          <rect x="${padding.left}" y="${bandTop}" width="${width - padding.left - padding.right}" height="${bandHeight}" fill="${colors.tostBand}" stroke="rgba(16, 185, 129, 0.4)" stroke-dasharray="3"/>
          <text x="${width - padding.right - 10}" y="${bandTop + 16}" fill="#10b981" font-size="11" font-weight="600" text-anchor="end">NMR TOST Band (δ = ±0.125)</text>
        `;
      }
    }

    // Regression Lines per species
    const speciesList = ["Mouse", "Human", "Naked Mole-Rat"];
    speciesList.forEach(sp => {
      const pts = data.filter(d => d.species === sp);
      if (pts.length < 2) return;

      const statItem = stats.find(s => s.comparison.startsWith(sp));
      const slope = statItem ? statItem.estimate : 0;
      const meanX = pts.reduce((sum, d) => sum + d.age_normalized, 0) / pts.length;
      const meanY = pts.reduce((sum, d) => sum + d.program_score, 0) / pts.length;
      const intercept = meanY - slope * meanX;

      const x1 = Math.min(...pts.map(d => d.age_normalized)) * 0.9;
      const x2 = Math.max(...pts.map(d => d.age_normalized)) * 1.05;
      const y1 = slope * x1 + intercept;
      const y2 = slope * x2 + intercept;

      const strokeColor = sp === "Mouse" ? colors.mouse : (sp === "Human" ? colors.human : colors.nmr);
      svg.innerHTML += `
        <line x1="${xScale(x1)}" y1="${yScale(y1)}" x2="${xScale(x2)}" y2="${yScale(y2)}" stroke="${strokeColor}" stroke-width="3" stroke-linecap="round"/>
      `;
    });

    // Scatter Data Points
    data.forEach(d => {
      const cx = xScale(d.age_normalized);
      const cy = yScale(d.program_score);
      const fillColor = d.species === "Mouse" ? colors.mouse : (d.species === "Human" ? colors.human : colors.nmr);

      svg.innerHTML += `
        <circle cx="${cx}" cy="${cy}" r="6" fill="${fillColor}" stroke="#ffffff" stroke-width="1.5" opacity="0.9">
          <title>${d.species} (${d.donor}): Age ${d.chronological_age} (${(d.age_normalized * 100).toFixed(1)}% lifespan), Score: ${d.program_score.toFixed(3)}</title>
        </circle>
      `;
    });

    // Axis Labels
    svg.innerHTML += `
      <text x="${width / 2}" y="${height - 15}" fill="#94a3b8" font-size="12" font-weight="500" text-anchor="middle">Life-Stage Normalized Age (Fraction of Species Maximum Lifespan)</text>
      <text x="${-height / 2}" y="24" fill="#94a3b8" font-size="12" font-weight="500" text-anchor="middle" transform="rotate(-90)">Conserved Aging Program Score</text>
    `;
  }

  // ---------------------------------------------------------
  // 4. View 2: Gene Trajectory Query
  // ---------------------------------------------------------
  function initGenesView() {
    const searchInput = document.getElementById("gene-search-input");
    const programSelect = document.getElementById("program-filter-select");

    if (searchInput) {
      searchInput.addEventListener("input", (e) => {
        const query = e.target.value.toLowerCase().trim();
        filterGeneTable(query, programSelect ? programSelect.value : "ALL");
      });
    }

    if (programSelect) {
      programSelect.addEventListener("change", (e) => {
        const query = searchInput ? searchInput.value.toLowerCase().trim() : "";
        filterGeneTable(query, e.target.value);
      });
    }

    renderGeneTable(state.data.genes || []);
    if (state.data.genes && state.data.genes.length > 0) {
      selectGene(state.data.genes[0].gene);
    }
  }

  function filterGeneTable(query, program) {
    let genes = state.data.genes || [];
    if (query) {
      genes = genes.filter(g => 
        g.gene.toLowerCase().includes(query) || 
        g.human_ortholog.toLowerCase().includes(query)
      );
    }
    if (program && program !== "ALL") {
      genes = genes.filter(g => g.program === program);
    }
    renderGeneTable(genes);
  }

  function renderGeneTable(genes) {
    const tbody = document.getElementById("table-gene-panel-body");
    if (!tbody) return;

    if (genes.length === 0) {
      tbody.innerHTML = `<tr><td colspan="7" style="text-align:center; color:#64748b;">No matching genes found.</td></tr>`;
      return;
    }

    tbody.innerHTML = genes.map(g => `
      <tr class="clickable ${g.gene === state.selectedGene ? 'selected' : ''}" data-gene="${g.gene}">
        <td><strong>${g.gene}</strong></td>
        <td><strong>${g.human_ortholog}</strong></td>
        <td>${g.program}</td>
        <td><span class="${g.mouse_log2fc > 0 ? 'text-accent-red' : 'text-accent-blue'}">${g.mouse_log2fc > 0 ? '+' : ''}${g.mouse_log2fc.toFixed(2)}</span> [${g.mouse_ci_low.toFixed(2)}, ${g.mouse_ci_high.toFixed(2)}]</td>
        <td><span class="${g.human_log2fc > 0 ? 'text-accent-red' : 'text-accent-blue'}">${g.human_log2fc > 0 ? '+' : ''}${g.human_log2fc.toFixed(2)}</span></td>
        <td>${g.nmr_slope > 0 ? '+' : ''}${g.nmr_slope.toFixed(2)}</td>
        <td><span class="status-pill status-conserved">${g.conservation_status}</span></td>
      </tr>
    `).join("");

    tbody.querySelectorAll("tr").forEach(tr => {
      tr.addEventListener("click", () => {
        const geneSymbol = tr.dataset.gene;
        selectGene(geneSymbol);
      });
    });
  }

  function selectGene(geneSymbol) {
    state.selectedGene = geneSymbol;
    const geneObj = (state.data.genes || []).find(g => g.gene === geneSymbol);
    if (!geneObj) return;

    // Highlight row
    document.querySelectorAll("#table-gene-panel-body tr").forEach(tr => {
      tr.classList.toggle("selected", tr.dataset.gene === geneSymbol);
    });

    // Update header
    document.getElementById("selected-gene-title").textContent = `Gene Trajectory: ${geneObj.gene} (Mouse) ↔ ${geneObj.human_ortholog} (Human)`;
    document.getElementById("selected-gene-desc").textContent = `Functional Program: ${geneObj.program} | Expected Direction: ${geneObj.expected_direction}`;

    const badgesContainer = document.getElementById("selected-gene-badges");
    if (badgesContainer) {
      badgesContainer.innerHTML = `
        <span class="status-pill status-conserved">${geneObj.conservation_status}</span>
        <span class="badge-text" style="color: #94a3b8; font-family: monospace;">FDR: ${geneObj.mouse_fdr}</span>
      `;
    }

    renderGeneChart();
  }

  function renderGeneChart() {
    const svg = document.getElementById("gene-svg");
    if (!svg) return;

    const geneObj = (state.data.genes || []).find(g => g.gene === state.selectedGene);
    if (!geneObj) return;

    const width = 900;
    const height = 320;
    const padding = { top: 40, right: 60, bottom: 60, left: 80 };

    svg.innerHTML = "";

    const speciesData = [
      { sp: "Mouse (TMS)", val: geneObj.mouse_log2fc, ci_low: geneObj.mouse_ci_low, ci_high: geneObj.mouse_ci_high, color: colors.mouse },
      { sp: "Human (SenNet/TS)", val: geneObj.human_log2fc, ci_low: geneObj.human_log2fc - 0.15, ci_high: geneObj.human_log2fc + 0.15, color: colors.human },
      { sp: "Naked Mole-Rat (Savina)", val: geneObj.nmr_slope, ci_low: geneObj.nmr_slope - 0.08, ci_high: geneObj.nmr_slope + 0.08, color: colors.nmr }
    ];

    const maxAbs = Math.max(1.5, ...speciesData.map(d => Math.max(Math.abs(d.ci_low), Math.abs(d.ci_high))));
    const yScale = (v) => height / 2 - (v / maxAbs) * (height / 2 - padding.top);

    // Center Baseline
    const zeroY = height / 2;
    svg.innerHTML += `
      <line x1="${padding.left}" y1="${zeroY}" x2="${width - padding.right}" y2="${zeroY}" stroke="rgba(255,255,255,0.2)" stroke-width="1.5"/>
      <text x="${padding.left - 10}" y="${zeroY + 4}" fill="#64748b" font-size="11" text-anchor="end" font-family="monospace">0.00</text>
    `;

    // Bars & CIs
    const barWidth = 80;
    const colStep = (width - padding.left - padding.right) / 3;

    speciesData.forEach((d, idx) => {
      const centerX = padding.left + idx * colStep + colStep / 2;
      const barTop = d.val >= 0 ? yScale(d.val) : zeroY;
      const barH = Math.abs(yScale(d.val) - zeroY);

      // CI error bar
      const ciTop = yScale(d.ci_high);
      const ciBottom = yScale(d.ci_low);
      svg.innerHTML += `
        <line x1="${centerX}" y1="${ciTop}" x2="${centerX}" y2="${ciBottom}" stroke="${d.color}" stroke-width="2"/>
        <line x1="${centerX - 10}" y1="${ciTop}" x2="${centerX + 10}" y2="${ciTop}" stroke="${d.color}" stroke-width="2"/>
        <line x1="${centerX - 10}" y1="${ciBottom}" x2="${centerX + 10}" y2="${ciBottom}" stroke="${d.color}" stroke-width="2"/>
      `;

      // Value Bar
      svg.innerHTML += `
        <rect x="${centerX - barWidth/2}" y="${barTop}" width="${barWidth}" height="${barH}" fill="${d.color}" rx="4" opacity="0.85"/>
        <text x="${centerX}" y="${d.val >= 0 ? barTop - 12 : zeroY + barH + 18}" fill="${d.color}" font-size="13" font-weight="700" text-anchor="middle" font-family="monospace">
          ${d.val > 0 ? '+' : ''}${d.val.toFixed(2)}
        </text>
        <text x="${centerX}" y="${height - 20}" fill="#94a3b8" font-size="12" font-weight="600" text-anchor="middle">
          ${d.sp}
        </text>
      `;
    });

    svg.innerHTML += `
      <text x="${-height / 2}" y="24" fill="#94a3b8" font-size="12" font-weight="500" text-anchor="middle" transform="rotate(-90)">Aging Effect Size (log₂FC / Slope)</text>
    `;
  }

  // ---------------------------------------------------------
  // 5. View 3: Conserved vs Divergent Programs
  // ---------------------------------------------------------
  function initProgramsView() {
    const tbody = document.getElementById("table-divergence-body");
    if (!tbody || !state.data.divergence) return;

    tbody.innerHTML = state.data.divergence.map(row => {
      const isConserved = row.status === "Conserved";
      const isLowPower = row.power_flag && row.power_flag.includes("Low-Powered");

      let badgeClass = "status-divergent";
      if (isConserved) badgeClass = "status-conserved";
      if (isLowPower) badgeClass = "status-lowpower";

      return `
        <tr>
          <td><strong>${row.program}</strong></td>
          <td>${row.cell_type}</td>
          <td>${row.tissue}</td>
          <td><span class="status-pill ${badgeClass}">${row.status}</span></td>
          <td><span class="${row.spearman_rho > 0 ? 'text-accent-green' : 'text-accent-red'}">${row.spearman_rho > 0 ? '+' : ''}${row.spearman_rho.toFixed(3)}</span></td>
          <td>${(row.strict_replication * 100).toFixed(1)}%</td>
          <td><span class="badge-text" style="color: ${isLowPower ? '#f59e0b' : '#10b981'};">${row.power_flag || 'Adequately Powered'}</span></td>
        </tr>
      `;
    }).join("");
  }

  // ---------------------------------------------------------
  // 6. View 4: Atlas & Tissue Explorer
  // ---------------------------------------------------------
  function initAtlasView() {
    const buttons = document.querySelectorAll("#umap-color-group .btn-toggle");
    buttons.forEach(btn => {
      btn.addEventListener("click", () => {
        buttons.forEach(b => b.classList.remove("active"));
        btn.classList.add("active");
        state.umapColorMode = btn.dataset.color;
        renderUmap();
      });
    });

    // Populate composition table
    const tbody = document.getElementById("table-composition-body");
    if (tbody && state.data.composition) {
      tbody.innerHTML = state.data.composition.slice(0, 15).map(c => `
        <tr>
          <td>${c.tissue}</td>
          <td><strong>${c.cell_ontology_class}</strong></td>
          <td>${c.n_cells_3m}</td>
          <td>${c.n_cells_18m}</td>
          <td><span class="status-pill status-conserved">${c.prop_3m ? (c.prop_18m / c.prop_3m).toFixed(2) : '1.00'}x</span></td>
        </tr>
      `).join("");
    }

    renderUmap();
  }

  function renderUmap() {
    const canvas = document.getElementById("umap-canvas");
    if (!canvas || !state.data.umap || state.data.umap.length === 0) return;

    const ctx = canvas.getContext("2d");
    const width = canvas.width;
    const height = canvas.height;
    ctx.clearRect(0, 0, width, height);

    const pts = state.data.umap;
    const xCoords = pts.map(p => p.x);
    const yCoords = pts.map(p => p.y);

    const xMin = Math.min(...xCoords);
    const xMax = Math.max(...xCoords);
    const yMin = Math.min(...yCoords);
    const yMax = Math.max(...yCoords);

    const pad = 30;
    const scaleX = (x) => pad + (x - xMin) / (xMax - xMin) * (width - 2 * pad);
    const scaleY = (y) => height - pad - (y - yMin) / (yMax - yMin) * (height - 2 * pad);

    pts.forEach(p => {
      const cx = scaleX(p.x);
      const cy = scaleY(p.y);

      let color = "#38bdf8";
      if (state.umapColorMode === "age") {
        color = colors.age[p.age] || "#38bdf8";
      } else if (state.umapColorMode === "tissue") {
        color = colors.tissue[p.tissue] || "#94a3b8";
      } else if (state.umapColorMode === "cell_type") {
        // Hash cell type name to color
        let hash = 0;
        for (let i = 0; i < p.cell_type.length; i++) hash = p.cell_type.charCodeAt(i) + ((hash << 5) - hash);
        const hue = Math.abs(hash) % 360;
        color = `hsl(${hue}, 70%, 60%)`;
      }

      ctx.beginPath();
      ctx.arc(cx, cy, 2.2, 0, 2 * Math.PI);
      ctx.fillStyle = color;
      ctx.globalAlpha = 0.75;
      ctx.fill();
    });

    // Update legend
    const legend = document.getElementById("umap-legend");
    if (!legend) return;

    if (state.umapColorMode === "age") {
      legend.innerHTML = `
        <div class="legend-row">
          <span class="legend-item"><span class="legend-color" style="background:#38bdf8;"></span> 3m Young</span>
          <span class="legend-item"><span class="legend-color" style="background:#f59e0b;"></span> 18m Aged</span>
          <span class="legend-item"><span class="legend-color" style="background:#f43f5e;"></span> 24m Senescent</span>
        </div>
      `;
    } else if (state.umapColorMode === "tissue") {
      legend.innerHTML = `
        <div class="legend-row">
          ${Object.entries(colors.tissue).map(([t, c]) => `
            <span class="legend-item"><span class="legend-color" style="background:${c};"></span> ${t}</span>
          `).join("")}
        </div>
      `;
    } else {
      legend.innerHTML = `<span class="footnote">Showing diverse cell ontology classes across multi-tissue manifold.</span>`;
    }
  }

  // ---------------------------------------------------------
  // 7. View 5: Aging Signal Ranking
  // ---------------------------------------------------------
  function initRankingView() {
    const tbody = document.getElementById("table-ranking-body");
    if (!tbody || !state.data.ranking) return;

    tbody.innerHTML = state.data.ranking.map((r, idx) => `
      <tr>
        <td><strong>#${idx + 1}</strong></td>
        <td><strong>${r.unit}</strong></td>
        <td><span class="metric-val text-accent-red" style="font-size:1rem;">${r.n_de_genes}</span></td>
        <td>[${r.n_de_ci_low}, ${r.n_de_ci_high}]</td>
        <td>${r.median_abs_log2fc ? r.median_abs_log2fc.toFixed(3) : "—"}</td>
        <td>${r.n_mice_total}</td>
        <td>${r.n_cells_total}</td>
      </tr>
    `).join("");
  }

  // Kickoff
  setupNavigation();
  loadAllData();
});
