/**
 * Agescape Webapp (Spec v2.0): Main Application Engine
 * Renders G0 through G13 via Vendored ECharts 5.5.1 + Direct DOM/SVG
 */

(function () {
  "use strict";

  // Application State
  const state = {
    data: {},
    selectedUnit: null,
    powerReadout: {
      donors: 4,
      sd: 8,
      dispersion: 0.1,
      min_30y_log2fc: 1.474,
      verdict: "LOW-POWERED"
    },
    activeCharts: {}
  };

  // Pinned Color System
  const PALETTE = {
    mouse: "#2563eb",
    human: "#7c3aed",
    nmr: "#ea580c",
    dirUp: "#dc2626",
    dirDown: "#475569",
    tissueLung: "#0284c7",
    tissueKidney: "#0d9488",
    tissueNeutral: "#94a3b8",
    tostBand: "rgba(234, 88, 12, 0.18)",
    tostStroke: "rgba(234, 88, 12, 0.6)",
    verdictAgree: "#0284c7",
    verdictPartial: "#f59e0b"
  };

  function createChart(dom) {
    const existing = echarts.getInstanceByDom(dom);
    if (existing) existing.dispose();
    const chart = echarts.init(dom);
    const setOption = chart.setOption.bind(chart);
    chart.setOption = (option, ...args) => {
      if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) option.animation = false;
      return setOption(option, ...args);
    };
    return chart;
  }

  // --------------------------------------------------------------------------
  // Initialization & Data Loading
  // --------------------------------------------------------------------------
  async function init() {
    setupModalsAndDrawers();
    setupUnitPanel();
    try {
      const files = [
        "atlas", "calibration", "agreement", "design-grid",
        "program", "signature", "savina", "power-grid", "anage", "humancohorts", "dlpfc"
      ];
      const fetches = files.map(f => fetch(`data/${f}.json`).then(r => {
        if (!r.ok) throw new Error(`HTTP ${r.status} on ${f}.json`);
        return r.json();
      }));

      const loaded = await Promise.all(fetches);
      files.forEach((f, idx) => {
        state.data[f] = loaded[idx];
      });

      console.log("Agescape app data loaded successfully:", Object.keys(state.data));

      setupNavigation();

      // Render Charts
      renderG0AtlasMap();
      renderG1RankingBars();
      renderG4MiniGrid();
      renderG3Agreement();
      renderG2Calibration();
      renderG4FullGrid();
      renderG6ProgramBars();
      renderG7Matrix();
      renderG9SignatureMultiples();
      renderG10HonestInterval();
      renderG11LifespanContext();
      renderG14AnAge();
      renderG15HumanCohorts();
      renderG16LifespanAtlas();
      renderG12SavinaCurve();
      renderG13PowerMap();
      setupAccessibilityEnhancements();

      // Check URL hash for initial unit panel open
      handleUrlHash();
      window.addEventListener("hashchange", handleUrlHash);

      let resizeTimer;
      let compactLayout = window.innerWidth < 560;
      window.addEventListener("resize", () => {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(() => {
          if ((window.innerWidth < 560) !== compactLayout) {
            compactLayout = window.innerWidth < 560;
            [renderG0AtlasMap, renderG1RankingBars, renderG3Agreement,
              renderG6ProgramBars, renderG9SignatureMultiples,
              renderG10HonestInterval, renderG15HumanCohorts,
              renderG16LifespanAtlas].forEach(render => {
                const domId = { renderG0AtlasMap: "g0", renderG1RankingBars: "g1" }[render.name];
                if (domId && state.activeCharts[domId]) state.activeCharts[domId].dispose();
                render();
              });
          }
          Object.values(state.activeCharts).forEach(c => c && c.resize && !c.isDisposed() && c.resize());
        }, 150);
      });

    } catch (err) {
      console.error("Failed to initialize Agescape webapp:", err);
      const notice = document.createElement("p");
      notice.className = "load-error";
      notice.setAttribute("role", "alert");
      notice.textContent = "The interactive charts could not load. Reload the page, or use the Code link to inspect the analysis.";
      document.querySelector(".hero-inner").appendChild(notice);
    }
  }

  // --------------------------------------------------------------------------
  // Navigation & Scroll Tracking
  // --------------------------------------------------------------------------
  function setupNavigation() {
    const chapters = document.querySelectorAll(".chapter");
    const railSteps = document.querySelectorAll(".rail-step");

    const observer = new IntersectionObserver((entries) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          const chId = entry.target.id;
          const num = chId.replace("ch", "");
          railSteps.forEach(s => {
            s.classList.toggle("active", s.dataset.chapter === num);
          });
        }
      });
    }, { rootMargin: "-20% 0px -60% 0px" });

    chapters.forEach(c => observer.observe(c));

    // Intro hero: activate the "Start" rail dot while the hero is in view
    const intro = document.getElementById("intro");
    if (intro) {
      const introObserver = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
          if (entry.isIntersecting) {
            railSteps.forEach(s => {
              s.classList.toggle("active", s.dataset.chapter === "0");
            });
          }
        });
      }, { rootMargin: "-20% 0px -60% 0px" });
      introObserver.observe(intro);
    }
  }

  // --------------------------------------------------------------------------
  // Modals & Drawers
  // --------------------------------------------------------------------------
  function setupModalsAndDrawers() {
    // Methods modal
    const methodsBtn = document.getElementById("btn-open-methods");
    const methodsModal = document.getElementById("methods-modal");
    const methodsClose = document.getElementById("btn-close-methods");

    const closeMethods = () => {
      methodsModal.hidden = true;
      document.body.style.overflow = "";
      document.querySelector("main").inert = false;
      document.querySelector(".top-nav").inert = false;
      methodsBtn.focus();
    };
    if (methodsBtn && methodsModal && methodsClose) {
      methodsBtn.addEventListener("click", () => {
        methodsModal.hidden = false;
        document.body.style.overflow = "hidden";
        document.querySelector("main").inert = true;
        document.querySelector(".top-nav").inert = true;
        methodsClose.focus();
        Object.values(state.activeCharts).forEach(c => c.dispatchAction && c.dispatchAction({ type: "hideTip" }));
      });
      methodsClose.addEventListener("click", closeMethods);
      methodsModal.addEventListener("click", (e) => {
        if (e.target === methodsModal) closeMethods();
      });
      methodsModal.addEventListener("keydown", (e) => {
        if (e.key === "Escape") { e.preventDefault(); closeMethods(); }
        if (e.key !== "Tab") return;
        const items = [...methodsModal.querySelectorAll("button, a[href], input, [tabindex='0']")];
        const first = items[0], last = items[items.length - 1];
        if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
        else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
      });
    }

    // Provenance Drawers
    document.querySelectorAll(".drawer-toggle").forEach(btn => {
      btn.addEventListener("click", () => {
        const dId = btn.dataset.drawer;
        const drawer = document.getElementById(dId);
        if (drawer) {
          const isHidden = drawer.hidden;
          drawer.hidden = !isHidden;
          btn.setAttribute("aria-expanded", String(isHidden));
        }
      });
    });

    document.querySelectorAll(".drawer-close").forEach(btn => {
      btn.addEventListener("click", () => {
        const drawer = btn.closest(".drawer");
        if (drawer) {
          drawer.hidden = true;
          const toggle = document.querySelector(`[data-drawer="${drawer.id}"]`);
          if (toggle) toggle.setAttribute("aria-expanded", "false");
        }
      });
    });

    // Compare Published toggles
    document.querySelectorAll(".btn-compare").forEach(btn => {
      btn.setAttribute("aria-controls", btn.dataset.target);
      btn.setAttribute("aria-expanded", "false");
      btn.addEventListener("click", () => {
        const comp = document.getElementById(btn.dataset.target);
        if (comp) {
          comp.hidden = !comp.hidden;
          btn.setAttribute("aria-expanded", String(!comp.hidden));
        }
      });
    });

    document.querySelectorAll(".compare-close").forEach(btn => {
      btn.setAttribute("aria-label", "Close published figure comparison");
      btn.addEventListener("click", () => {
        const comp = btn.closest(".compare-container");
        if (comp) {
          comp.hidden = true;
          const toggle = document.querySelector(`[data-target="${comp.id}"]`);
          if (toggle) { toggle.setAttribute("aria-expanded", "false"); toggle.focus(); }
        }
      });
    });

    // CSV Download handlers
    document.querySelectorAll(".btn-download").forEach(btn => {
      btn.setAttribute("aria-label", btn.textContent.trim() + " (CSV)");
      btn.addEventListener("click", () => {
        handleCsvDownload(btn.dataset.export, btn.dataset.sub);
      });
    });
  }

  // --------------------------------------------------------------------------
  // CSV Download Generator
  // --------------------------------------------------------------------------
  function handleCsvDownload(exportName, subKey) {
    const rawData = state.data[exportName];
    if (!rawData) return;

    let rows = [];
    let filename = `${exportName}.csv`;

    if (exportName === "atlas") {
      rows = subKey === "ranked" ? rawData.ranked : rawData.units;
      filename = `atlas-${subKey || "units"}.csv`;
    } else if (exportName === "calibration") {
      rows = subKey === "posctl" ? rawData.positive_controls : [{
        unit: "Lung::adventitial cell", observed: rawData.observed,
        permutation_mean: rawData.perm_mean, permutation_replicates: 5,
        within_3m_sex_de_genes: rawData.sex_floor
      }];
      filename = subKey === "posctl" ? "calibration-positive-controls.csv" : "calibration-summary.csv";
    } else if (exportName === "agreement") {
      if (subKey === "scatter") {
        rows = rawData.scatter.points;
        filename = `agreement-kidney-pct-scatter.csv`;
      } else {
        rows = rawData.lollipop;
        filename = `agreement-lollipop-units.csv`;
      }
    } else if (exportName === "program") {
      rows = subKey === "catalog" ? rawData.units : rawData.genes;
      filename = `program-${subKey || "genes"}.csv`;
    } else if (exportName === "signature") {
      if (subKey === "points") {
        rows = Object.entries(rawData.arms).flatMap(([arm, value]) =>
          value.points.map(point => ({ arm, ...point })));
        filename = `signature-all-donor-scores.csv`;
      } else {
        rows = [
          { test: "mouse_minus_nmr_slope_difference", ...rawData.stats.diff },
          { test: "nmr_tost_equivalence", ...rawData.stats.tost }
        ];
        filename = `signature-statistical-summary.csv`;
      }
    } else if (exportName === "savina") {
      rows = rawData.curve;
      filename = `savina-reproduction-curve.csv`;
    } else if (exportName === "power-grid") {
      rows = rawData.tiles;
      filename = `human-power-simulation-grid.csv`;
    } else if (exportName === "anage") {
      rows = Array.isArray(rawData) ? rawData : [];
      filename = `anage-species-longevity.csv`;
    } else if (exportName === "humancohorts") {
      rows = subKey === "celltypes" ? (rawData.celltypes || []) : (rawData.cohorts || []);
      filename = subKey === "celltypes" ? `human-cohort-celltype-regressions.csv` : `human-cohort-summary.csv`;
    } else if (exportName === "dlpfc") {
      rows = rawData.celltypes || [];
      filename = `dlpfc-lifespan-atlas-celltype-regressions.csv`;
    } else if (exportName === "design-grid") {
      const dg = rawData;
      rows = dg.age_groups.map((ag, i) => ({
        age_group: ag,
        cells_female: dg.cells.female[i],
        cells_male: dg.cells.male[i],
        mice_female: dg.mice.female[i],
        mice_male: dg.mice.male[i]
      }));
      filename = `sex-age-design-grid.csv`;
    }

    if (!rows || rows.length === 0) return;

    downloadCsv(rows, filename);
  }

  function downloadCsv(rows, filename) {
    const headers = [...new Set(rows.flatMap(row => Object.keys(row)))];
    const cell = value => {
      const text = value === undefined || value === null ? "" : String(value);
      return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
    };
    const csvContent = [headers.map(cell).join(","),
      ...rows.map(row => headers.map(h => cell(row[h])).join(","))].join("\n");
    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  // --------------------------------------------------------------------------
  // Accessibility: chart text twins & Escape-to-close
  // --------------------------------------------------------------------------
  function setupAccessibilityEnhancements() {
    document.querySelectorAll(".echart-container, .mini-chart").forEach(dom => {
      const block = dom.closest(".graph-block") || dom.closest(".expander-content") || dom.closest("section");
      const titleEl = block ? block.querySelector(".graph-title, .sub-chart-title") : null;
      const label = titleEl ? titleEl.textContent.trim() : "Aging data chart";
      dom.setAttribute("role", "img");
      dom.setAttribute("aria-label", `Chart: ${label}. Use the adjacent description and CSV download for the values.`);
    });

    const unitSelect = document.createElement("select");
    unitSelect.className = "accessible-chart-select";
    unitSelect.setAttribute("aria-label", "Inspect a testable mouse cell-type unit");
    unitSelect.innerHTML = "<option value=''>Inspect a cell-type unit</option>" +
      state.data.atlas.ranked.map(r => `<option value="${r.unit}">${r.unit.replace("::", ": ")} (${r.n_de_genes} genes)</option>`).join("");
    unitSelect.addEventListener("change", () => { if (unitSelect.value) openUnitPanel(unitSelect.value); });
    document.getElementById("block-g1").appendChild(unitSelect);

    const powerSelect = document.createElement("select");
    powerSelect.className = "accessible-chart-select";
    powerSelect.setAttribute("aria-label", "Select a study design and update the power readout");
    powerSelect.innerHTML = "<option value=''>Select a study design</option>" +
      state.data["power-grid"].tiles.map((t, i) => `<option value="${i}">${t.n_donors} donors, age SD ${t.sd}, dispersion ${t.dispersion}: effect ${t.min_30y_log2fc}</option>`).join("");
    powerSelect.addEventListener("change", () => {
      if (powerSelect.value !== "") updatePowerReadout(state.data["power-grid"].tiles[Number(powerSelect.value)]);
    });
    document.getElementById("power-readout").before(powerSelect);

    document.addEventListener("keydown", (e) => {
      if (e.key !== "Escape") return;
      const panel = document.getElementById("unit-panel");
      if (panel && panel.classList.contains("open")) { closeUnitPanel(); return; }
      const openDrawer = document.querySelector(".drawer:not([hidden])");
      if (openDrawer) {
        openDrawer.hidden = true;
        const toggle = document.querySelector(`[data-drawer="${openDrawer.id}"]`);
        if (toggle) { toggle.setAttribute("aria-expanded", "false"); toggle.focus(); }
        return;
      }
      const openCompare = document.querySelector(".compare-container:not([hidden])");
      if (openCompare) {
        openCompare.hidden = true;
        const toggle = document.querySelector(`[data-target="${openCompare.id}"]`);
        if (toggle) { toggle.setAttribute("aria-expanded", "false"); toggle.focus(); }
      }
    });
  }

  // --------------------------------------------------------------------------
  // Slide-Over Unit Panel (Section 8)
  // --------------------------------------------------------------------------
  function setupUnitPanel() {
    const panel = document.getElementById("unit-panel");
    const closeBtn = document.getElementById("panel-close-btn");
    const searchInput = document.getElementById("panel-gene-search");
    const downloadBtn = document.getElementById("btn-download-panel-genes");
    panel.inert = true;
    panel.setAttribute("role", "dialog");
    panel.setAttribute("aria-modal", "true");
    panel.setAttribute("aria-labelledby", "panel-unit-name");

    document.getElementById("panel-human-link").addEventListener("click", closeUnitPanel);
    panel.addEventListener("keydown", (e) => {
      if (e.key !== "Tab") return;
      const items = [...panel.querySelectorAll("button, a[href], input, select")];
      const first = items[0], last = items[items.length - 1];
      if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
      else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
    });

    if (closeBtn) {
      closeBtn.addEventListener("click", () => {
        closeUnitPanel();
      });
    }

    if (searchInput) {
      searchInput.addEventListener("input", (e) => {
        filterUnitGenes(e.target.value.toLowerCase().trim());
      });
    }

    if (downloadBtn) {
      downloadBtn.addEventListener("click", () => {
        if (!state.selectedUnit) return;
        const prog = state.data.program;
        if (!prog) return;
        const genes = prog.genes.filter(g => g.unit === state.selectedUnit);
        downloadCsv(genes, `${state.selectedUnit.replace("::", "-")}-genes.csv`);
      });
    }
  }

  function openUnitPanel(unitName) {
    const panel = document.getElementById("unit-panel");
    if (!panel || (state.selectedUnit === unitName && panel.classList.contains("open"))) return;
    state.panelReturnFocus = document.activeElement;
    state.selectedUnit = unitName;

    panel.inert = false;
    panel.classList.add("open");
    panel.setAttribute("aria-hidden", "false");
    document.getElementById("panel-gene-search").value = "";

    const pcb = document.getElementById("panel-close-btn");
    if (pcb) pcb.focus();

    const parts = unitName.split("::");
    const tissue = parts[0];
    const ct = parts[1] || "";

    document.getElementById("panel-tissue").textContent = tissue;
    document.getElementById("panel-unit-name").textContent = ct;

    // Unit program stats
    const prog = state.data.program;
    const unitCat = prog ? prog.units.find(u => u.unit === unitName) : null;
    const nUp = unitCat ? unitCat.n_up : 0;
    const nDown = unitCat ? unitCat.n_down : 0;
    const nTotal = nUp + nDown;

    document.getElementById("panel-stat-up").textContent = nUp;
    document.getElementById("panel-stat-down").textContent = nDown;
    document.getElementById("panel-stat-total").textContent = nTotal;

    // Agreement & guards
    const agr = state.data.agreement;
    const agrItem = agr ? agr.lollipop.find(u => u.unit === unitName) : null;
    if (agrItem) {
      document.getElementById("panel-meta-agreement").textContent = `ρ = ${agrItem.spearman_rho.toFixed(3)} (${agrItem.verdict})`;
    } else {
      document.getElementById("panel-meta-agreement").textContent = "N/A";
    }

    if (unitCat) {
      document.getElementById("panel-meta-housekeeping").textContent = `${unitCat.frac_ribo_mito.toFixed(3)} (background ${unitCat.frac_ribo_mito_background.toFixed(3)})`;
    }

    // Render Mini Stacked Bar in Panel
    renderPanelBar(nUp, nDown);

    // Populate Gene Table
    renderUnitGeneTable(unitName);

    // Update URL Hash
    const slug = unitName.toLowerCase().replace(/[^a-z0-9]+/g, "-");
    window.location.hash = `unit-${slug}`;
  }

  function closeUnitPanel() {
    const panel = document.getElementById("unit-panel");
    if (panel) {
      panel.classList.remove("open");
      panel.setAttribute("aria-hidden", "true");
      panel.inert = true;
    }
    state.selectedUnit = null;
    if (window.location.hash.startsWith("#unit-")) {
      history.replaceState("", document.title, window.location.pathname + window.location.search);
    }
    if (state.panelReturnFocus && state.panelReturnFocus.isConnected) state.panelReturnFocus.focus();
  }

  function handleUrlHash() {
    const hash = window.location.hash;
    if (hash.startsWith("#unit-")) {
      const slug = hash.replace("#unit-", "");
      const atlas = state.data.atlas;
      if (atlas) {
        const found = atlas.ranked.find(r => r.unit.toLowerCase().replace(/[^a-z0-9]+/g, "-") === slug);
        if (found) openUnitPanel(found.unit);
      }
    } else if (state.selectedUnit) {
      closeUnitPanel();
    }
  }

  function renderPanelBar(nUp, nDown) {
    const dom = document.getElementById("chart-panel-bar");
    if (!dom) return;

    let chart = state.activeCharts["panel-bar"];
    if (!chart) {
      chart = createChart(dom);
      state.activeCharts["panel-bar"] = chart;
    }

    chart.setOption({
      grid: { left: 10, right: 10, top: 15, bottom: 20 },
      xAxis: { type: "value", show: false },
      yAxis: { type: "category", data: [""], show: false },
      series: [
        {
          name: "Upregulated",
          type: "bar",
          stack: "total",
          data: [nUp],
          itemStyle: { color: PALETTE.dirUp, borderRadius: [4, 0, 0, 4] }
        },
        {
          name: "Downregulated",
          type: "bar",
          stack: "total",
          data: [nDown],
          itemStyle: { color: PALETTE.dirDown, borderRadius: [0, 4, 4, 0] }
        }
      ]
    });
  }

  function renderUnitGeneTable(unitName) {
    const prog = state.data.program;
    if (!prog) return;

    const genes = prog.genes.filter(g => g.unit === unitName);
    const tbody = document.getElementById("panel-gene-tbody");
    if (!tbody) return;

    tbody.innerHTML = genes.map(g => `
      <tr>
        <td><strong>${g.gene}</strong></td>
        <td><span style="color:${g.direction === 'up' ? PALETTE.dirUp : '#64748b'}">${g.direction.toUpperCase()}</span></td>
        <td class="font-mono">${g.log2fc !== null ? g.log2fc.toFixed(2) : "N/A"}</td>
        <td class="font-mono">${g.fdr !== null ? g.fdr.toFixed(4) : "N/A"}</td>
      </tr>
    `).join("");
  }

  function filterUnitGenes(query) {
    const rows = document.querySelectorAll("#panel-gene-tbody tr");
    rows.forEach(tr => {
      const geneName = tr.cells[0].textContent.toLowerCase();
      tr.style.display = geneName.includes(query) ? "" : "none";
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 1: G0 Atlas Eligibility Map
  // --------------------------------------------------------------------------
  function renderG0AtlasMap() {
    const dom = document.getElementById("chart-g0");
    const atlas = state.data.atlas;
    if (!dom || !atlas) return;

    const chart = createChart(dom);
    state.activeCharts["g0"] = chart;

    // Group units by tissue
    const tissueMap = {};
    atlas.units.forEach(u => {
      if (!tissueMap[u.tissue]) tissueMap[u.tissue] = [];
      tissueMap[u.tissue].push(u);
    });

    const tissues = Object.keys(tissueMap).sort((a, b) => {
      const aPow = tissueMap[a].filter(u => u.powered).length;
      const bPow = tissueMap[b].filter(u => u.powered).length;
      return bPow - aPow; // Lung and Kidney first
    });

    // Custom scatter layout representing tiles
    const data = [];
    tissues.forEach((tissue, rowIdx) => {
      const unitsInTissue = tissueMap[tissue];
      unitsInTissue.forEach((u, colIdx) => {
        data.push({
          value: [colIdx, rowIdx, u.powered ? 1 : 0],
          unit: u.unit,
          tissue: u.tissue,
          ct: u.cell_ontology_class,
          powered: u.powered
        });
      });
    });

    const humanize = s => s.replace(/_/g, " ").replace(/\b\w/g, c => c.toUpperCase()).replace(/ And /g, " and ");
    const yAxisLabels = tissues.map(t => {
      const powCount = tissueMap[t].filter(u => u.powered).length;
      const totCount = tissueMap[t].length;
      return `${humanize(t)} (${powCount}/${totCount})`;
    });

    chart.setOption({
      tooltip: {
        trigger: "item",
        formatter: (params) => {
          const d = params.data;
          const status = d.powered ? "POWERED (Eligible)" : "Ineligible (&lt;3 mice or &lt;50 cells)";
          return `<strong>${d.tissue}</strong><br>${d.ct}<br><span style="color:${d.powered ? '#0369a1' : '#64748b'}">${status}</span>`;
        }
      },
      grid: dom.clientWidth < 560
        ? { left: 104, right: 10, top: 10, bottom: 20 }
        : { left: 160, right: 30, top: 10, bottom: 20 },
      xAxis: {
        type: "value",
        show: false,
        min: -0.5,
        max: Math.max(...Object.values(tissueMap).map(units => units.length)) - 0.5
      },
      yAxis: {
        type: "category",
        data: yAxisLabels,
        inverse: true,
        axisLine: { lineStyle: { color: "rgba(15,23,42,0.2)" } },
        axisLabel: { color: "#64748b", fontSize: dom.clientWidth < 560 ? 9 : 11,
          width: dom.clientWidth < 560 ? 94 : null, overflow: "truncate" }
      },
      series: [{
        type: "scatter",
        symbol: "rect",
        symbolSize: dom.clientWidth < 560 ? [5, 7] : [11, 11],
        data: data,
        itemStyle: {
          color: (p) => {
            if (!p.data.powered) return "transparent";
            return p.data.tissue === "Lung" ? PALETTE.tissueLung : PALETTE.tissueKidney;
          },
          borderColor: (p) => {
            if (p.data.powered) {
              return p.data.tissue === "Lung" ? PALETTE.tissueLung : PALETTE.tissueKidney;
            }
            return "rgba(15,23,42,0.3)";
          },
          borderWidth: 1.5
        }
      }]
    });

    chart.on("click", (params) => {
      if (params.data && params.data.powered) {
        openUnitPanel(params.data.unit);
      }
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 1: G1 Ranking Bars
  // --------------------------------------------------------------------------
  function renderG1RankingBars() {
    const dom = document.getElementById("chart-g1");
    const atlas = state.data.atlas;
    if (!dom || !atlas) return;

    const chart = createChart(dom);
    state.activeCharts["g1"] = chart;

    const ranked = [...atlas.ranked].reverse(); // Ascending for horizontal bar
    const isMobile = dom.clientWidth < 560;
    const categories = ranked.map(r => {
      const full = r.unit.replace("::", " · ");
      if (isMobile) return full.length > 38 ? full.slice(0, 37) + "…" : full;
      const sep = full.indexOf(" · ");
      return sep > 0 ? full.slice(0, sep) + " ·\n" + full.slice(sep + 3) : full;
    });
    const values = ranked.map(r => r.n_de_genes);

    chart.setOption({
      tooltip: {
        trigger: "axis",
        axisPointer: { type: "shadow" },
        formatter: (params) => {
          const idx = params[0].dataIndex;
          const r = ranked[idx];
          return `<strong>${r.unit.replace("::", " · ")}</strong><br>Aging genes: <strong>${r.n_de_genes}</strong><br>Median |log₂FC|: ${r.median_abs_log2fc.toFixed(3)}`;
        }
      },
      grid: isMobile
        ? { left: 110, right: 46, top: 15, bottom: 25 }
        : { left: 300, right: 60, top: 15, bottom: 25 },
      xAxis: {
        type: "value",
        axisLabel: { color: "#64748b", fontSize: 11 },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: categories,
        axisLabel: { color: "#64748b", fontSize: 11 }
      },
      series: [{
        type: "bar",
        data: values,
        label: {
          show: true,
          position: "right",
          color: "#475569",
          fontFamily: "monospace",
          fontSize: 11
        },
        itemStyle: {
          color: (p) => {
            const tissue = ranked[p.dataIndex].tissue;
            return tissue === "Lung" ? PALETTE.tissueLung : PALETTE.tissueKidney;
          },
          borderRadius: [0, 4, 4, 0]
        }
      }]
    });

    chart.on("click", (params) => {
      const idx = params.dataIndex;
      openUnitPanel(ranked[idx].unit);
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 1 & 2: G4 Design Grid (Mini & Full)
  // --------------------------------------------------------------------------
  function renderG4MiniGrid() {
    renderG4GridImpl("chart-g4-mini", [40, 20, 10, 20]);
  }

  function renderG4FullGrid() {
    renderG4GridImpl("chart-g4-full", [60, 40, 20, 30]);
  }

  function renderG4GridImpl(domId, gridPadding) {
    const dom = document.getElementById(domId);
    const dg = state.data["design-grid"];
    if (!dom || !dg) return;

    const chart = createChart(dom);
    state.activeCharts[domId] = chart;

    chart.setOption({
      tooltip: {
        trigger: "item",
        formatter: (params) => {
          const d = params.data;
          return `Age: <strong>${dg.age_groups[d[0]]}</strong><br>Sex: <strong>${d[1] === 0 ? 'Female' : 'Male'}</strong><br>Mice: <strong>${d[2]}</strong>`;
        }
      },
      grid: { left: gridPadding[0], right: gridPadding[1], top: gridPadding[2], bottom: gridPadding[3] },
      xAxis: {
        type: "category",
        data: dg.age_groups,
        axisLabel: { color: "#64748b" }
      },
      yAxis: {
        type: "category",
        data: ["Female", "Male"],
        axisLabel: { color: "#64748b" }
      },
      series: [{
        type: "heatmap",
        data: [
          ...dg.mice.female.map((cnt, i) => [i, 0, cnt]),
          ...dg.mice.male.map((cnt, i) => [i, 1, cnt])
        ],
        label: {
          show: true,
          color: "#ffffff",
          fontFamily: "monospace",
          fontWeight: "bold",
          formatter: (p) => p.data[2] > 0 ? `${p.data[2]} mice` : "0"
        },
        itemStyle: {
          color: (p) => {
            const cnt = p.data[2];
            const isFairFight = p.data[0] === 1 || p.data[0] === 2; // 3m and 18m
            if (cnt === 0) return "rgba(15,23,42,0.04)";
            return isFairFight ? "#0284c7" : "rgba(2, 132, 199, 0.35)";
          },
          borderWidth: 2,
          borderColor: (p) => (p.data[0] === 1 || p.data[0] === 2) ? "#0369a1" : "rgba(15,23,42,0.25)"
        }
      }]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 2: G3 Agreement (Scatter + Lollipop)
  // --------------------------------------------------------------------------
  function renderG3Agreement() {
    renderG3Scatter();
    renderG3Lollipop();
  }

  function renderG3Scatter() {
    const dom = document.getElementById("chart-g3-scatter");
    const agr = state.data.agreement;
    if (!dom || !agr || !agr.scatter) return;

    const chart = createChart(dom);
    state.activeCharts["g3-scatter"] = chart;

    const sc = agr.scatter;
    const pts = sc.points.map(p => [p.ours_bridge_slope, p.paper_coef, p.gene]);

    chart.setOption({
      tooltip: {
        trigger: "item",
        formatter: (params) => {
          const d = params.data;
          return `<strong>${d[2]}</strong><br>Our Bridge Slope (per month): ${d[0].toFixed(4)}<br>Paper Age Coef: ${d[1].toFixed(4)}`;
        }
      },
      title: {
        text: `Kidney proximal convoluted tubule: Spearman rho = ${sc.rho.toFixed(3)} across ${sc.n.toLocaleString()} genes`,
        left: 60,
        top: 10,
        textStyle: { color: "#0369a1", fontSize: 13, fontWeight: 600 }
      },
      grid: { left: 60, right: 30, top: 45, bottom: 40 },
      xAxis: {
        type: "value",
        name: "Our Bridge Slope",
        nameLocation: "middle",
        nameGap: 24,
        nameTextStyle: { color: "#64748b" },
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "value",
        name: "TMS Paper Age Coef",
        nameLocation: "middle",
        nameGap: 36,
        nameTextStyle: { color: "#64748b" },
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      series: [{
        type: "scatter",
        symbolSize: 3,
        data: pts,
        itemStyle: { color: "#0369a1", opacity: 0.4 }
      }]
    });
  }

  function renderG3Lollipop() {
    const dom = document.getElementById("chart-g3-lollipop");
    const agr = state.data.agreement;
    if (!dom || !agr) return;

    const chart = createChart(dom);
    state.activeCharts["g3-lollipop"] = chart;

    const units = [...agr.lollipop].reverse();
    const isMobile3 = dom.clientWidth < 560;
    const categories = units.map(u => {
      const full = u.unit.replace("::", " · ");
      if (isMobile3) return full.length > 30 ? full.slice(0, 29) + "…" : full;
      const sep = full.indexOf(" · ");
      return sep > 0 ? full.slice(0, sep) + " ·\n" + full.slice(sep + 3) : full;
    });
    const values = units.map(u => u.spearman_rho);

    chart.setOption({
      tooltip: {
        trigger: "axis",
        formatter: (p) => {
          const idx = p[0].dataIndex;
          const u = units[idx];
          const conc = u.sign_concordance != null
            ? `<br>Sign concordance (paper's significant genes): ${(u.sign_concordance * 100).toFixed(1)}%` : "";
          return `<strong>${u.unit}</strong><br>Spearman ρ: ${u.spearman_rho.toFixed(3)}${conc}<br>Verdict: <strong>${u.verdict}</strong>`;
        }
      },
      grid: isMobile3
        ? { left: 120, right: 30, top: 10, bottom: 25 }
        : { left: 300, right: 40, top: 10, bottom: 25 },
      xAxis: {
        type: "value",
        min: 0,
        max: 1.0,
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: categories,
        axisLabel: { color: "#64748b", fontSize: 10, lineHeight: 12 }
      },
      series: [{
        type: "scatter",
        symbolSize: 8,
        data: values,
        itemStyle: {
          color: (p) => units[p.dataIndex].verdict === "agree" ? PALETTE.verdictAgree : PALETTE.verdictPartial
        }
      }]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 2: G2 Calibration Dot-Plot
  // --------------------------------------------------------------------------
  function renderG2Calibration() {
    const dom = document.getElementById("chart-g2");
    const cal = state.data.calibration;
    if (!dom || !cal) return;

    const chart = createChart(dom);
    state.activeCharts["g2"] = chart;

    chart.setOption({
      grid: { left: 40, right: 60, top: 40, bottom: 40 },
      xAxis: {
        type: "value",
        min: -2,
        max: 60,
        axisLabel: { color: "#64748b", fontFamily: "monospace" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: ["DE Genes"],
        show: false
      },
      series: [
        {
          type: "scatter",
          symbolSize: 14,
          data: [
            { value: [cal.sex_floor, 0],
              name: "Differential expression by sex, 3-month group (2 female vs 5 male mice):\n0 genes at FDR < 0.05",
              itemStyle: { color: "rgba(0,0,0,0)", borderColor: "#64748b", borderWidth: 2 },
              label: { show: true, position: "top", distance: 6, align: "left", formatter: "Differential expression by sex, 3-month group (2 female vs 5 male mice):\n0 genes at FDR < 0.05", color: "#475569", fontSize: 11 } },
            { value: [cal.perm_mean, 0],
              name: "Label-permutation null, mean: 4.6",
              itemStyle: { color: "#f59e0b" },
              label: { show: true, position: "top", distance: 22, align: "left", formatter: "Label-permutation null, mean: 4.6", color: "#475569", fontSize: 11 } },
            { value: [cal.observed, 0],
              name: "Observed, 3 month vs 18 month: 55",
              itemStyle: { color: "#0369a1" },
              label: { show: true, position: "top", distance: 6, align: "right", formatter: "Observed, 3 month vs 18 month: 55", color: "#475569", fontSize: 11 } }
          ]
        }
      ]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 3: G6 Program Direction Bars
  // --------------------------------------------------------------------------
  function renderG6ProgramBars() {
    const dom = document.getElementById("chart-g6");
    const prog = state.data.program;
    if (!dom || !prog) return;

    const chart = createChart(dom);
    state.activeCharts["g6"] = chart;

    const units = [...prog.units].reverse();
    const isMobile6 = dom.clientWidth < 560;
    const categories = units.map(u => {
      const full = u.unit.replace("::", " · ");
      if (isMobile6) return full.length > 30 ? full.slice(0, 29) + "…" : full;
      const sep = full.indexOf(" · ");
      return sep > 0 ? full.slice(0, sep) + " ·\n" + full.slice(sep + 3) : full;
    });

    chart.setOption({
      tooltip: {
        trigger: "axis",
        axisPointer: { type: "shadow" },
        formatter: (p) => {
          const idx = p[0].dataIndex;
          const u = units[idx];
          return `<strong>${u.unit}</strong><br>Up with age: <strong style="color:${PALETTE.dirUp}">${u.n_up}</strong><br>Down with age: <strong style="color:${PALETTE.dirDown}">${u.n_down}</strong><br>Total: ${u.n_up + u.n_down}`;
        }
      },
      grid: isMobile6
        ? { left: 120, right: 40, top: 15, bottom: 25 }
        : { left: 300, right: 60, top: 15, bottom: 25 },
      xAxis: {
        type: "value",
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: categories,
        axisLabel: { color: "#64748b", fontSize: 11, lineHeight: 13 }
      },
      series: [
        {
          name: "Up",
          type: "bar",
          stack: "total",
          data: units.map(u => u.n_up),
          itemStyle: { color: PALETTE.dirUp }
        },
        {
          name: "Down",
          type: "bar",
          stack: "total",
          data: units.map(u => u.n_down),
          itemStyle: { color: PALETTE.dirDown, borderRadius: [0, 4, 4, 0] }
        }
      ]
    });

    chart.on("click", (p) => {
      openUnitPanel(units[p.dataIndex].unit);
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 3: G7 Gene × Unit Dot Matrix
  // --------------------------------------------------------------------------
  function renderG7Matrix() {
    const dom = document.getElementById("chart-g7");
    const prog = state.data.program;
    if (!dom || !prog || !prog.matrix_top40) return;

    const chart = createChart(dom);
    state.activeCharts["g7"] = chart;

    const m40 = prog.matrix_top40;
    const geneNames = m40.rows.map(r => r.gene);
    const unitNames = m40.columns.map(u => u.replace("::", " · "));

    // Map entries to scatter points
    const data = [];
    prog.genes.forEach(g => {
      const gIdx = geneNames.indexOf(g.gene);
      const uIdx = m40.columns.indexOf(g.unit);
      if (gIdx >= 0 && uIdx >= 0) {
        data.push([uIdx, gIdx, g.direction === "up" ? 1 : -1, g.gene, g.unit, g.log2fc]);
      }
    });

    chart.setOption({
      tooltip: {
        trigger: "item",
        formatter: (p) => {
          const d = p.data;
          return `Gene: <strong>${d[3]}</strong><br>Unit: ${d[4]}<br>Direction: <strong>${d[2] === 1 ? 'UP' : 'DOWN'}</strong> (log₂FC = ${d[5].toFixed(2)})`;
        }
      },
      grid: { left: 90, right: 30, top: 40, bottom: 140 },
      xAxis: {
        type: "category",
        data: unitNames,
        axisLabel: { color: "#64748b", rotate: 35, fontSize: 10 }
      },
      yAxis: {
        type: "category",
        data: geneNames,
        inverse: true,
        axisLabel: { color: "#64748b", fontSize: 10, fontFamily: "monospace" }
      },
      series: [{
        type: "scatter",
        symbolSize: 8,
        data: data,
        itemStyle: {
          color: (p) => p.data[2] === 1 ? PALETTE.dirUp : PALETTE.dirDown
        }
      }]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 4: G9 Signature Small Multiples
  // --------------------------------------------------------------------------
  function renderG9SignatureMultiples() {
    const dom = document.getElementById("chart-g9");
    const sig = state.data.signature;
    if (!dom || !sig) return;

    const chart = createChart(dom);
    state.activeCharts["g9"] = chart;

    const arms = sig.arms;
    const mousePts = arms.mouse_wb.points.map(p => [p.age_normalized, p.score, p.donor]);
    const humanPts = arms.human_skin ? arms.human_skin.points.map(p => [p.age_normalized, p.score, p.donor]) : [];
    const nmrPts = arms.nmr ? arms.nmr.points.map(p => [p.age_normalized, p.score, p.donor]) : [];
    const isMobile = dom.clientWidth < 560;

    // Trend line through the committed points, mean-centered over the
    // observed age range: y = mean(y) + slope * (x - mean(x)).
    const fitLine = (pts, slope) => {
      if (!pts.length) return [];
      const mx = pts.reduce((a, p) => a + p[0], 0) / pts.length;
      const my = pts.reduce((a, p) => a + p[1], 0) / pts.length;
      const xs = pts.map(p => p[0]);
      const lo = Math.min(...xs), hi = Math.max(...xs);
      return [[lo, my + slope * (lo - mx)], [hi, my + slope * (hi - mx)]];
    };

    chart.setOption({
      title: [
        { text: isMobile ? `Mouse (n = ${arms.mouse_wb.n})` : `Mouse (slope = ${arms.mouse_wb.slope.toFixed(3)}, p = ${arms.mouse_wb.p.toFixed(3)}, n = ${arms.mouse_wb.n})`, left: "3%", top: "2%", textStyle: { color: PALETTE.mouse, fontSize: 12 } },
        { text: isMobile ? `Human (n = ${arms.human_skin.n}, low-powered)` : `Human (slope = ${arms.human_skin.slope.toFixed(3)}, n = ${arms.human_skin.n}): LOW-POWERED`, left: "3%", top: "35%", textStyle: { color: PALETTE.human, fontSize: 12 } },
        { text: isMobile ? `Naked mole-rat (n = ${arms.nmr.n})` : `Naked Mole-Rat (slope = ${arms.nmr.slope.toFixed(3)}, p = ${arms.nmr.p.toFixed(3)}, n = ${arms.nmr.n})`, left: "3%", top: "68%", textStyle: { color: PALETTE.nmr, fontSize: 12 } }
      ],
      grid: [
        { left: 50, right: 30, top: "7%", height: "23%" },
        { left: 50, right: 30, top: "40%", height: "23%" },
        { left: 50, right: 30, top: "73%", height: "23%" }
      ],
      xAxis: [
        { gridIndex: 0, min: 0, max: 1.0, show: false },
        { gridIndex: 1, min: 0, max: 1.0, show: false },
        { gridIndex: 2, min: 0, max: 1.0, axisLabel: { color: "#64748b" }, name: "Normalized Lifespan (0 to 1)", nameLocation: "middle", nameGap: 24, nameTextStyle: { color: "#64748b" } }
      ],
      yAxis: [
        { gridIndex: 0, min: -0.5, max: 1.5, axisLabel: { color: "#64748b" } },
        { gridIndex: 1, min: -0.5, max: 1.5, axisLabel: { color: "#64748b" } },
        { gridIndex: 2, min: -0.5, max: 1.5, axisLabel: { color: "#64748b" } }
      ],
      series: [
        // Mouse
        {
          xAxisIndex: 0, yAxisIndex: 0,
          type: "scatter",
          symbolSize: 7,
          data: mousePts,
          itemStyle: { color: PALETTE.mouse }
        },
        {
          xAxisIndex: 0, yAxisIndex: 0,
          type: "line",
          data: fitLine(mousePts, arms.mouse_wb.slope),
          itemStyle: { color: PALETTE.mouse },
          lineStyle: { width: 2.5 }
        },
        // Human
        {
          xAxisIndex: 1, yAxisIndex: 1,
          type: "scatter",
          symbolSize: 7,
          data: humanPts,
          itemStyle: { color: PALETTE.human }
        },
        {
          xAxisIndex: 1, yAxisIndex: 1,
          type: "line",
          data: fitLine(humanPts, arms.human_skin.slope),
          itemStyle: { color: PALETTE.human },
          lineStyle: { type: "dashed", width: 2 }
        },
        // NMR
        {
          xAxisIndex: 2, yAxisIndex: 2,
          type: "scatter",
          symbolSize: 7,
          data: nmrPts,
          itemStyle: { color: PALETTE.nmr }
        },
        {
          xAxisIndex: 2, yAxisIndex: 2,
          type: "line",
          data: fitLine(nmrPts, arms.nmr.slope),
          itemStyle: { color: PALETTE.nmr },
          lineStyle: { width: 2.5 }
        }
      ]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 4: G15 Human cohort arm (M9) — forest plot of program slopes
  // --------------------------------------------------------------------------
  function renderG15HumanCohorts() {
    const dom = document.getElementById("chart-g15");
    const hc = state.data.humancohorts;
    if (!dom || !hc) return;

    const chart = createChart(dom);
    state.activeCharts["g15"] = chart;

    const HUMAN = "#7c3aed"; // single human accent colour; hollow/filled carries the verdict
    const humanize = n => n.replace(/_/g, " ").replace("ImmuneAging", "Immune-Aging");
    const rows = hc.cohorts.map((c, i) => ({
      name: `${humanize(c.cohort)} (${c.n_donors})`,
      short: humanize(c.cohort),
      value: [c.slope_per_year * 1000, i],
      ci: [c.ci_lo * 1000, c.ci_hi * 1000],
      n: c.n_donors, p: c.p, q: c.q_bh, group: c.group, verdict: c.verdict,
      significant: c.ci_hi < 0 || c.ci_lo > 0,
    }));

    chart.setOption({
      tooltip: {
        trigger: "item",
        formatter: (p) => {
          const d = p.data;
          return `<strong>${d.short}</strong> (${d.group})<br>donors: ${d.n}, slope: ${d.value[0].toFixed(3)}x10^-3/yr<br>95% CI: [${d.ci[0].toFixed(3)}, ${d.ci[1].toFixed(3)}]<br>p = ${d.p.toExponential(2)}, BH q = ${d.q.toExponential(2)}${d.significant ? " (nominal CI excludes 0)" : ""}`;
        }
      },
      grid: dom.clientWidth < 560
        ? { left: 112, right: 14, top: 30, bottom: 55 }
        : { left: 185, right: 60, top: 30, bottom: 40 },
      xAxis: {
        type: "value",
        name: dom.clientWidth < 560 ? "score slope (x10^-3/year)" : "program score slope per year (x10^-3)",
        nameLocation: "middle",
        nameGap: 28,
        min: -3.5, max: 3.5, interval: 1,
        axisLabel: { color: "#64748b", fontSize: dom.clientWidth < 560 ? 9 : 11 },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: rows.map(r => r.name),
        axisLabel: { color: "#334155", fontSize: dom.clientWidth < 560 ? 9 : 11,
          width: dom.clientWidth < 560 ? 102 : null, overflow: "truncate" },
        axisTick: { show: false }, axisLine: { show: false }
      },
      series: [
        ...rows.map(d => ({
          name: `${d.short} 95% CI`,
          type: "line",
          data: [[d.ci[0], d.value[1]], [d.ci[1], d.value[1]]],
          symbol: "none",
          silent: true,
          lineStyle: { color: HUMAN, width: d.significant ? 3 : 2, opacity: d.significant ? 1 : 0.55 },
          tooltip: { show: false }
        })),
        {
          name: "Cohort estimate",
          type: "scatter",
          data: rows.map(d => ({
            ...d,
            itemStyle: { color: d.significant ? HUMAN : "#ffffff", borderColor: HUMAN,
              borderWidth: 2, opacity: d.significant ? 1 : 0.85 }
          })),
          symbolSize: (value, params) => params.data.significant ? 13 : 11,
          label: {
            show: true,
            position: "top",
            distance: 7,
            formatter: d => d.data.significant ? `p=${d.data.p.toFixed(3)}; q=${d.data.q.toFixed(3)}` : "",
            color: "#0f172a", fontSize: 10, fontWeight: 700
          },
          markLine: {
            silent: true, symbol: "none",
            lineStyle: { color: "#0f172a", width: 2, type: "solid" },
            data: [{ xAxis: 0 }],
            label: { show: true, formatter: "no change", position: "end", color: "#0f172a", fontWeight: 600 }
          }
        }
      ]
    });

    // top cell-type hits table (exploratory ranking, not a confirmed hit list)
    const tbody = document.getElementById("g15-tbody");
    if (tbody) {
      const fmt = v => (v * 1000).toFixed(2);
      const NEURON_DAGGER = (r) => (r.cohort === "MSSM_DLPFC" && /L2\/3 intratelencephalic/.test(r.cell_type))
        ? " &#8224;" : "";
      tbody.innerHTML = hc.celltypes.slice(0, 10).map(r => `
        <tr>
          <td>${humanize(r.cohort)}</td>
          <td>${r.cell_type}${NEURON_DAGGER(r)}</td>
          <td class="font-mono">${r.n_donors}</td>
          <td class="font-mono">${r.slope_per_year >= 0 ? "+" : ""}${fmt(r.slope_per_year)}</td>
          <td class="font-mono">${r.p.toExponential(2)}</td>
          <td>${r.verdict === "POWERED-FLOOR"
            ? `<span class="badge-sc">meets floor</span>`
            : `<span class="species-none">low power</span>`}</td>
        </tr>`).join("");
    }
  }

  // --------------------------------------------------------------------------
  // CHAPTER 5: G16 DLPFC Lifespan Atlas (M10 human lifespan arm)
  // --------------------------------------------------------------------------
  function renderG16LifespanAtlas() {
    const dom = document.getElementById("chart-g16");
    const dl = state.data.dlpfc;
    if (!dom || !dl) return;

    const chart = createChart(dom);
    state.activeCharts["g16"] = chart;

    const PURPLE = "#7c3aed";       // human cohort colour family (matches G15 brain)
    const PURPLE_SOFT = "#6d28d9";  // CPM robustness arm, AA-visible on white
    const x = v => v * 1000;        // slopes shown x10^-3/yr
    const isMobile = dom.clientWidth < 560;

    // sort by primary slope descending; glia rise, immune fall
    const rows = [...dl.celltypes].sort((a, b) => b.slope_per_year - a.slope_per_year)
      .map((c, i) => ({ ...c, row: i }));

    chart.setOption({
      title: {
        text: "", left: "center", top: 2,
        subtext: isMobile
          ? `available-cell-type average: p = ${dl.cohort.whole_cortex.p.toFixed(2)}`
          : `available-cell-type average: ${x(dl.cohort.whole_cortex.slope_per_year).toFixed(2)}x10^-3/yr, p = ${dl.cohort.whole_cortex.p.toFixed(2)}; no rise detected`,
        subtextStyle: { fontSize: 12, fontWeight: 700, color: "#0f172a" }
      },
      tooltip: {
        trigger: "item",
        formatter: (p) => {
          const d = rows[p.dataIndex];
          const verdict = d.concordant ? "same-direction nominal signal in both schemes"
            : d.direction === "normalization-dependent direction" ? "direction reverses between schemes: not claimed"
            : d.direction === "scheme-sensitive significance" ? "same sign, but significance is not robust across schemes"
            : "no nominal signal in either scheme";
          return `<strong>${d.cell_type}</strong> (${d.n_donors} donors)<br>` +
            `raw-sum: ${x(d.slope_per_year).toFixed(2)}x10^-3/yr, p = ${d.p.toExponential(2)}, BH q = ${d.q_bh.toExponential(2)}<br>` +
            `CPM: ${x(d.slope_cpm).toFixed(2)}x10^-3/yr, p = ${d.p_cpm.toExponential(2)}, BH q = ${d.q_bh_cpm.toExponential(2)}<br>` +
            `adults-only p = ${d.p_adult.toExponential(2)}<br>${verdict}`;
        }
      },
      grid: isMobile
        ? { left: 92, right: 16, top: 42, bottom: 40 }
        : { left: 118, right: 110, top: 42, bottom: 40 },
      xAxis: {
        type: "value",
        name: "program score slope per year (x10^-3)",
        nameLocation: "middle",
        nameGap: 28,
        min: -3.5, max: 3.5, interval: 1,
        axisLabel: { color: "#64748b", fontSize: isMobile ? 10 : 11 },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: rows.map(r => r.cell_type),
        inverse: true,
        axisLabel: { color: "#334155", fontSize: isMobile ? 10 : 11 },
        axisTick: { show: false }, axisLine: { show: false }
      },
      series: [
        {
          type: "custom",
          data: rows,
          renderItem: (params, api) => {
            const d = rows[params.dataIndex];
            const y = api.coord([0, d.row])[1];
            const xLo = api.coord([x(d.ci_lo), d.row])[0];
            const xHi = api.coord([x(d.ci_hi), d.row])[0];
            const xM = api.coord([x(d.slope_per_year), d.row])[0];
            const xC = api.coord([x(d.slope_cpm), d.row])[0];
            const children = [];
            // raw-sum CI whisker + dot (primary)
            children.push({ type: "line", shape: { x1: xLo, y1: y, x2: xHi, y2: y },
              style: { stroke: PURPLE, lineWidth: d.concordant ? 3 : 2, opacity: d.concordant ? 1 : 0.5 } });
            children.push(d.concordant
              ? { type: "circle", shape: { cx: xM, cy: y, r: 6 },
                  style: { fill: PURPLE, stroke: "#0f172a", lineWidth: 1 } }
              : { type: "circle", shape: { cx: xM, cy: y, r: 5 },
                  style: { fill: "rgba(0,0,0,0)", stroke: PURPLE, lineWidth: 2, opacity: 0.65 } });
            // CPM robustness dot, smaller, nudged below the row center
            children.push(d.concordant
              ? { type: "circle", shape: { cx: xC, cy: y + 9, r: 4 }, style: { fill: PURPLE_SOFT } }
              : { type: "circle", shape: { cx: xC, cy: y + 9, r: 3.5 },
                  style: { fill: "rgba(0,0,0,0)", stroke: PURPLE_SOFT, lineWidth: 1.6, opacity: 0.7 } });
            // direct direction label for concordant rows (desktop only; on mobile the tooltip carries it)
            if (d.concordant && !isMobile) {
              const label = d.direction === "up" ? "rises" : "falls";
              const gridLeft = isMobile ? 92 : 118;
              let lx, align;
              if (d.slope_per_year >= 0) { lx = xHi + 8; align = "left"; }
              else if (xLo - 46 > gridLeft) { lx = xLo - 8; align = "right"; }
              else { lx = xM + 12; align = "left"; }
              children.push({ type: "text", style: {
                x: lx, y: y - 6,
                text: label, textAlign: align,
                fontSize: 11, fontWeight: 700, fill: "#0f172a" } });
            }
            return { type: "group", children };
          }
        },
        {
          type: "scatter", data: [],
          markLine: {
            silent: true, symbol: "none",
            lineStyle: { color: "#0f172a", width: 2 },
            data: [{ xAxis: 0, label: { show: true, formatter: "no change", position: "end", color: "#0f172a", fontWeight: 600 } },
              { xAxis: x(dl.cohort.whole_cortex.slope_per_year),
                lineStyle: { color: "#94a3b8", width: 1.5, type: "dashed" },
                label: { show: false } }]
          }
        }
      ]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 4: G10 The Honest Interval (Forest Plot)
  // --------------------------------------------------------------------------
  function renderG10HonestInterval() {
    const dom = document.getElementById("chart-g10");
    const sig = state.data.signature;
    if (!dom || !sig) return;

    const chart = createChart(dom);
    state.activeCharts["g10"] = chart;

    const delta = sig.stats.tost.delta;

    chart.setOption({
      grid: dom.clientWidth < 560
        ? { left: 82, right: 16, top: 25, bottom: 35 }
        : { left: 160, right: 60, top: 25, bottom: 35 },
      xAxis: {
        type: "value",
        min: -0.4,
        max: 0.6,
        axisLabel: { color: "#64748b", fontFamily: "monospace" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: ["Slope Difference", "NMR Slope", "Mouse Slope"],
        axisLabel: { color: "#64748b", fontSize: dom.clientWidth < 560 ? 9 : 11,
          width: dom.clientWidth < 560 ? 72 : null, overflow: "truncate" }
      },
      series: [
        // Bold line at 0
        {
          type: "line",
          markLine: {
            silent: true,
            symbol: "none",
            data: [{ xAxis: 0, lineStyle: { color: "rgba(15,23,42,0.6)", width: 2 } }]
          }
        },
        // TOST equivalence band (NMR row) + CI whiskers + direct label (custom layer)
        {
          type: "custom",
          z: 1,
          renderItem: (params, api) => {
            const kind = api.value(4);
            if (kind === 2) { // band
              const p1 = api.coord([api.value(0), api.value(2)]);
              const p2 = api.coord([api.value(1), api.value(2)]);
              const h = 16;
              return {
                type: "rect",
                shape: { x: p1[0], y: p1[1] - h, width: p2[0] - p1[0], height: 2 * h },
                style: { fill: PALETTE.tostBand, stroke: PALETTE.tostStroke, lineWidth: 1 }
              };
            }
            if (kind === 3) { // direct label
              const p = api.coord([api.value(0), api.value(2)]);
              const shortLabel = dom.clientWidth < 480;
              return {
                type: "text",
                style: {
                  x: p[0], y: p[1], text: shortLabel ? "equivalence unproven" : "equivalence not established",
                  fill: "#9a3412", fontSize: shortLabel ? 9 : 11, fontWeight: 600,
                  align: "left", verticalAlign: "middle"
                }
              };
            }
            // kind 1: CI whisker with end caps (row color from closure, not data dims)
            const p1 = api.coord([api.value(0), api.value(2)]);
            const p2 = api.coord([api.value(1), api.value(2)]);
            const y = p1[1];
            const c = api.value(2) === 2 ? PALETTE.mouse : (api.value(2) === 1 ? PALETTE.nmr : "#0369a1");
            return {
              type: "group",
              children: [
                { type: "line", shape: { x1: p1[0], y1: y, x2: p2[0], y2: y }, style: { stroke: c, lineWidth: 2 } },
                { type: "line", shape: { x1: p1[0], y1: y - 6, x2: p1[0], y2: y + 6 }, style: { stroke: c, lineWidth: 2 } },
                { type: "line", shape: { x1: p2[0], y1: y - 6, x2: p2[0], y2: y + 6 }, style: { stroke: c, lineWidth: 2 } }
              ]
            };
          },
          data: [
            { value: [-delta, delta, 1, 0, 2] },
            { value: [sig.stats.diff.ci[0], sig.stats.diff.ci[1], 0, 0, 1] },
            { value: [sig.arms.nmr.ci[0], sig.arms.nmr.ci[1], 1, 0, 1] },
            { value: [sig.arms.mouse_wb.ci[0], sig.arms.mouse_wb.ci[1], 2, 0, 1] },
            { value: [delta + 0.015, 1, 1, 0, 3] }
          ],
          tooltip: { show: false },
          silent: true
        },
        // Point estimates
        {
          type: "scatter",
          symbolSize: 10,
          z: 3,
          data: [
            [sig.stats.diff.slope, 0],
            [sig.arms.nmr.slope, 1],
            [sig.arms.mouse_wb.slope, 2]
          ],
          itemStyle: {
            color: (p) => p.dataIndex === 2 ? PALETTE.mouse : (p.dataIndex === 1 ? PALETTE.nmr : "#0369a1")
          }
        }
      ]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 4: G11 Lifespan Context
  // --------------------------------------------------------------------------
  function renderG11LifespanContext() {
    const dom = document.getElementById("chart-g11");
    const sig = state.data.signature;
    if (!dom || !sig || !sig.lifespans) return;

    const chart = createChart(dom);
    state.activeCharts["g11"] = chart;

    chart.setOption({
      grid: { left: 120, right: 60, top: 15, bottom: 25 },
      xAxis: {
        type: "value",
        max: 130,
        axisLabel: { color: "#64748b", formatter: "{value} y" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "category",
        data: ["Human", "Naked Mole-Rat", "Mouse"],
        axisLabel: { color: "#64748b" }
      },
      series: [
        {
          type: "bar",
          data: [
            { value: sig.lifespans.human, itemStyle: { color: PALETTE.human } },
            { value: sig.lifespans.nmr, itemStyle: { color: PALETTE.nmr } },
            { value: sig.lifespans.mouse, itemStyle: { color: PALETTE.mouse } }
          ],
          barWidth: 16,
          markArea: {
            silent: true,
            itemStyle: { color: "rgba(15,23,42,0.48)" },
            label: { show: true, color: "#ffffff", fontSize: 9, formatter: "tested" },
            data: [
              [{ xAxis: 1 / 12, yAxis: "Mouse" }, { xAxis: 2.5, yAxis: "Mouse" }],
              [{ xAxis: 1.1, yAxis: "Naked Mole-Rat" }, { xAxis: 11.3, yAxis: "Naked Mole-Rat" }]
            ]
          },
          label: {
            show: true,
            position: "right",
            color: "#64748b",
            formatter: (p) => `${p.value}y max`
          }
        }
      ]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 4: G14 AnAge Species Browser
  // --------------------------------------------------------------------------
  const SC_SPECIES = {
    "Homo sapiens": { color: "#7c3aed", label: "Human" },
    "Heterocephalus glaber": { color: "#ea580c", label: "Naked mole-rat" },
    "Mus musculus": { color: "#2563eb", label: "Mouse" }
  };

  function renderG14AnAge() {
    const dom = document.getElementById("chart-g14");
    const data = state.data.anage;
    if (!dom || !Array.isArray(data)) return;

    const chart = createChart(dom);
    state.activeCharts["g14"] = chart;

    const scatterPts = data.filter(d => d.kg && d.long).map(d => ({
      value: [d.kg, d.long],
      name: d.n,
      cls: d.cls,
      ref: d,
      itemStyle: d.sc === "yes" ? { color: SC_SPECIES[d.latin].color, borderColor: "#0f172a", borderWidth: 1.5 } : undefined
    }));

    const mammals = scatterPts.filter(p => p.cls === "Mammalia" && !p.itemStyle);
    const others = scatterPts.filter(p => p.cls !== "Mammalia" && !p.itemStyle);
    const featured = scatterPts.filter(p => p.itemStyle);

    chart.setOption({
      tooltip: {
        trigger: "item",
        formatter: (p) => {
          const d = p.data.ref;
          const latin = d ? `<br><em>${d.latin}</em>` : "";
          const mass = p.value[0] < 0.01 ? `${(p.value[0] * 1000).toFixed(0)} g` : `${p.value[0]} kg`;
          return `<strong>${p.data.name}</strong>${latin}<br>Class: ${p.data.cls}<br>Max longevity: <strong>${p.value[1]} yrs</strong><br>Adult mass: ${mass}`;
        }
      },
      grid: { left: 85, right: 30, top: 30, bottom: 55 },
      xAxis: {
        type: "log",
        name: "Adult body mass (kg, log scale)",
        nameLocation: "middle",
        nameGap: 30,
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "log",
        name: "Max longevity (yrs, log)",
        nameGap: 15,
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      series: [
        {
          type: "scatter",
          name: "Other classes",
          data: others,
          symbolSize: 5,
          itemStyle: { color: "rgba(100,116,139,0.35)" },
          emphasis: { scale: 1.5 }
        },
        {
          type: "scatter",
          name: "Mammals",
          data: mammals,
          symbolSize: 6,
          itemStyle: { color: "rgba(3,105,161,0.55)" },
          emphasis: { scale: 1.5 }
        },
        {
          type: "scatter",
          name: "Single-cell data available",
          data: featured,
          symbolSize: 15,
          label: {
            show: true,
            position: "top",
            color: "#0f172a",
            fontWeight: 700,
            formatter: (p) => {
              const d = p.data.ref;
              return (d && SC_SPECIES[d.latin] && SC_SPECIES[d.latin].label) || p.data.name;
            }
          },
          z: 10
        }
      ]
    });

    // Table with search
    const tbody = document.getElementById("species-tbody");
    const search = document.getElementById("species-search");
    const note = document.getElementById("species-count-note");
    if (!tbody) return;

    const sorted = [...data].sort((a, b) => b.long - a.long);
    const PAGE = 50;

    function esc(s) { return String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }

    function renderRows(filter) {
      const q = (filter || "").trim().toLowerCase();
      const hits = q
        ? sorted.filter(d => d.n.toLowerCase().includes(q) || d.latin.toLowerCase().includes(q))
        : sorted;
      const shown = hits.slice(0, PAGE);
      tbody.innerHTML = shown.map(d => `
        <tr>
          <td><strong>${esc(d.n)}</strong><br><span class="species-latin">${esc(d.latin)}</span></td>
          <td>${esc(d.cls)}</td>
          <td class="font-mono">${d.long}</td>
          <td class="font-mono">${d.kg != null ? d.kg.toLocaleString("en-US") : "n/a"}</td>
          <td>${d.sc === "yes"
            ? `<span class="badge-sc">In this explorer</span>`
            : `<span class="species-none">Not included</span>`}</td>
        </tr>`).join("");
      if (note) {
        note.textContent = q
          ? `${hits.length} matching species. Showing the first ${shown.length} by longevity.`
          : `Showing ${shown.length} of ${hits.length} species, sorted by maximum longevity. Search to filter; the CSV download holds the full table.`;
      }
    }

    renderRows("");
    if (search) {
      search.addEventListener("input", () => renderRows(search.value));
    }
  }

  // --------------------------------------------------------------------------
  // CHAPTER 4: G12 Savina Curve
  // --------------------------------------------------------------------------
  function renderG12SavinaCurve() {
    const dom = document.getElementById("chart-g12");
    const sav = state.data.savina;
    if (!dom || !sav) return;

    const chart = createChart(dom);
    state.activeCharts["g12"] = chart;

    const pts = sav.curve.map(r => [r.log2fc_cutoff, r.n_de_pval05]);

    chart.setOption({
      tooltip: {
        trigger: "axis",
        formatter: (p) => `Cutoff |log₂FC|: ${p[0].data[0]}<br>DE Genes: <strong>${p[0].data[1]}</strong>`
      },
      grid: { left: 60, right: 30, top: 25, bottom: 35 },
      xAxis: {
        type: "value",
        name: "log₂FC Cutoff",
        nameLocation: "middle",
        nameGap: 22,
        nameTextStyle: { color: "#64748b" },
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      yAxis: {
        type: "log",
        name: "DE Genes (log)",
        nameLocation: "middle",
        nameGap: 36,
        nameTextStyle: { color: "#64748b" },
        axisLabel: { color: "#64748b" },
        splitLine: { lineStyle: { color: "rgba(15,23,42,0.06)" } }
      },
      series: [
        {
          type: "line",
          data: pts,
          itemStyle: { color: PALETTE.nmr },
          lineStyle: { width: 2.5 }
        },
        {
          type: "scatter",
          symbolSize: 8,
          data: [
            { value: [1.0, 29], name: "Paper Rule: 29 genes" },
            { value: [1.45, 9], name: "Converted Rule: 9 genes (6 at FDR)" },
            { value: [0.25, 811], name: "Loosest cutoff: 811 by p rule (560 by FDR)" }
          ],
          itemStyle: { color: "#ffffff", borderColor: PALETTE.nmr, borderWidth: 2 }
        }
      ]
    });
  }

  // --------------------------------------------------------------------------
  // CHAPTER 5: G13 Power Map (3 Small Multiples)
  // --------------------------------------------------------------------------
  function renderG13PowerMap() {
    renderG13SinglePanel("chart-g13-low", 0.1);
    renderG13SinglePanel("chart-g13-mid", 0.2);
    renderG13SinglePanel("chart-g13-high", 0.4);
    updatePowerReadout(state.powerReadout);
  }

  function renderG13SinglePanel(domId, dispVal) {
    const dom = document.getElementById(domId);
    const pw = state.data["power-grid"];
    if (!dom || !pw) return;

    const chart = createChart(dom);
    state.activeCharts[domId] = chart;

    const tiles = pw.tiles.filter(t => Math.abs(t.dispersion - dispVal) < 0.01);
    const donorCats = pw.donors.map(String);
    const sdCats = pw.sds.map(String);

    // Color on a log scale: minimum detectable effect is multiplicative
    // (2.0 means twice as large an effect needed). One shared transform and
    // ramp across all three panels, so colours compare honestly panel to panel.
    const RAMP = ["#fef9c3", "#fde047", "#fbbf24", "#fb923c", "#ea580c", "#b91c1c", "#7f1d1d"];
    const hex = h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
    const rampColor = t => {
      const tt = Math.min(Math.max(t, 0), 1) * (RAMP.length - 1);
      const i = Math.min(Math.floor(tt), RAMP.length - 2);
      const f = tt - i;
      const a = hex(RAMP[i]), b = hex(RAMP[i + 1]);
      return `rgb(${a.map((v, k) => Math.round(v + (b[k] - v) * f)).join(",")})`;
    };
    const norm = v => (Math.log10(v) + 0.85) / 1.35; // 0.14 -> 0, 3.16 -> 1 (covers committed 0.167-2.949)

    const data = tiles.map(t => {
      const xIdx = pw.donors.indexOf(t.n_donors);
      const yIdx = pw.sds.indexOf(t.sd);
      return {
        value: [xIdx, yIdx, t.min_30y_log2fc, t.power_status, t.eligible_cohort, t.n_donors, t.sd, t.dispersion],
        itemStyle: { color: rampColor(norm(t.min_30y_log2fc)) },
        label: { color: t.min_30y_log2fc >= 1.35 ? "#ffffff" : "#3f3f46" }
      };
    });

    chart.setOption({
      tooltip: {
        trigger: "item",
        formatter: (p) => {
          const d = p.value;
          return `Donors: <strong>${d[5]}</strong> | SD: <strong>${d[6]}y</strong><br>Min 30y log₂FC: <strong>${d[2].toFixed(3)}</strong><br>Status: <strong>${d[4] ? "MEETS FLOOR" : "LOW-POWERED"}</strong>`;
        }
      },
      grid: { left: 40, right: 15, top: 15, bottom: 35 },
      xAxis: {
        type: "category",
        data: donorCats,
        name: "donors",
        nameLocation: "middle",
        nameGap: 22,
        nameTextStyle: { color: "#64748b", fontSize: 10 },
        axisLabel: { color: "#64748b", fontSize: 10 }
      },
      yAxis: {
        type: "category",
        data: sdCats,
        name: "age SD (years)",
        nameTextStyle: { color: "#64748b", fontSize: 10, align: "left" },
        axisLabel: { color: "#64748b", fontSize: 10 }
      },
      series: [{
        type: "heatmap",
        data: data,
        label: {
          show: true,
          fontSize: 9,
          fontFamily: "monospace",
          formatter: (p) => p.value[2].toFixed(2).replace(/^0/, "")
        },
        itemStyle: {
          borderWidth: 1.5,
          borderColor: (p) => p.value[4] ? "#0369a1" : "rgba(0,0,0,0.4)" // Outlined if meets floor
        },
        markPoint: Math.abs(dispVal - 0.1) < 0.01 ? {
          data: [{
            coord: [0, 0],
            value: "✕",
            itemStyle: { color: "rgba(0,0,0,0)" },
            label: { show: true, formatter: "✕", color: "#b91c1c", fontWeight: "bold", fontSize: 14, textBorderColor: "#ffffff", textBorderWidth: 2 }
          }]
        } : undefined
      }]
    });

    chart.on("click", (params) => {
      const d = params.value;
      updatePowerReadout({
        donors: d[5],
        sd: d[6],
        dispersion: d[7],
        min_30y_log2fc: d[2],
        verdict: d[3]
      });
    });
  }

  function updatePowerReadout(r) {
    const donors = r.donors ?? r.n_donors;
    const eligible = r.eligible_cohort ?? (r.verdict === "MEETS FLOOR" || r.verdict === "Adequately Powered");
    const verdict = eligible ? "MEETS FLOOR" : "LOW-POWERED";
    state.powerReadout = { ...r, donors, verdict };
    document.getElementById("readout-donors").textContent = donors;
    document.getElementById("readout-sd").textContent = `${r.sd} y`;
    document.getElementById("readout-disp").textContent = r.dispersion;
    document.getElementById("readout-log2fc").textContent = r.min_30y_log2fc.toFixed(3);

    const verdictEl = document.getElementById("readout-verdict");
    verdictEl.textContent = verdict;
    verdictEl.className = eligible ? "badge-verdict badge-floor" : "badge-verdict badge-lowpower";
  }

  // Kickoff on DOM Load
  document.addEventListener("DOMContentLoaded", init);
})();
