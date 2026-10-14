/* Conserved Aging Explorer — static, dependency-free vanilla JS.
 * Every displayed number comes from app/data/*.json, which notebooks/07-app-assets.ipynb
 * copies verbatim from results/*.csv (see data/manifest.json for source + SHA-256 per asset).
 * No external network calls: only same-origin fetches of the local JSON assets. */
(function () {
  "use strict";

  var DATA = {};
  var MANIFEST = null;

  var ASSET_NAMES = [
    "signature", "savina", "program_genes", "ranking",
    "gene_programs", "agreement", "confound", "umap_subsample"
  ];

  /* ---------------- formatting helpers ---------------- */
  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }
  function f(x, d) {
    if (x === null || x === undefined || Number.isNaN(x)) return "—";
    return Number(x).toFixed(d === undefined ? 3 : d);
  }
  function pfmt(p) {
    if (p === null || p === undefined || Number.isNaN(p)) return "—";
    p = Number(p);
    return p < 0.001 ? p.toExponential(2) : p.toFixed(3);
  }
  function ci(row) {
    if (row.ci_low === null || row.ci_low === undefined) return "—";
    return "[" + f(row.ci_low, 3) + ", " + f(row.ci_high, 3) + "]";
  }
  function byKey(rows, key) {
    var m = {};
    rows.forEach(function (r) { m[r[key]] = r; });
    return m;
  }

  /* ---------------- provenance ---------------- */
  function prov(assetKey) {
    var meta = MANIFEST && MANIFEST.assets && MANIFEST.assets[assetKey];
    if (!meta) return "";
    var srcs = meta.sources.map(function (s) { return "<code>" + esc(s.path) + "</code>"; }).join(" + ");
    return '<p class="provenance">Provenance: ' + srcs +
      ' — copied without modification by <code>notebooks/07-app-assets.ipynb</code>; ' +
      "SHA-256 of each source verified by <code>scripts/check_app_assets.py</code>.</p>";
  }

  /* ---------------- data loading ---------------- */
  function loadJSON(url) {
    return fetch(url).then(function (r) {
      if (!r.ok) throw new Error(url + " -> " + r.status);
      return r.json();
    });
  }
  function loadAll() {
    var jobs = { manifest: loadJSON("data/manifest.json") };
    ASSET_NAMES.forEach(function (n) { jobs[n] = loadJSON("data/" + n + ".json"); });
    var keys = Object.keys(jobs);
    return Promise.all(keys.map(function (k) { return jobs[k]; })).then(function (vals) {
      keys.forEach(function (k, i) {
        if (k === "manifest") MANIFEST = vals[i]; else DATA[k] = vals[i];
      });
    });
  }

  /* ================================================================
   * VIEW 1 — SIGNATURE (longevity-outlier contrast)
   * ================================================================ */
  function statsTableHtml() {
    var rows = DATA.signature.statistics.map(function (r) {
      var note = "";
      if (r.comparison.indexOf("TOST") !== -1) {
        note = r.p > 0.05
          ? ' <span class="chip bad">equivalence NOT established</span>'
          : ' <span class="chip ok">equivalence established</span>';
      }
      var hl = r.comparison.indexOf("TOST") !== -1 ? ' class="hl"' : "";
      return "<tr" + hl + "><td>" + esc(r.comparison) + note + "</td>" +
        '<td class="num">' + f(r.slope, 4) + "</td>" +
        '<td class="num">' + f(r.se, 4) + "</td>" +
        '<td class="num">' + ci(r) + "</td>" +
        '<td class="num">' + pfmt(r.p) + "</td>" +
        '<td class="num">' + (r.n == null ? "—" : r.n) + "</td></tr>";
    }).join("");
    return '<div class="table-wrap"><table class="data"><thead><tr>' +
      "<th>Comparison</th><th class='num'>Estimate</th><th class='num'>Std. err.</th>" +
      "<th class='num'>95% CI</th><th class='num'>p</th><th class='num'>n</th>" +
      "</tr></thead><tbody>" + rows + "</tbody></table></div>";
  }

  function gatesHtml() {
    var g = byKey(DATA.signature.gates, "gate");
    function cov(name) {
      var v = g[name] && g[name].value;
      return v == null ? "—" : (v * 100).toFixed(1) + "%";
    }
    var pc = g.mouse_positive_control;
    var cards = [
      { k: "Positive-control gate (mouse whole-body slope, 95% CI)", v: f(pc.value, 3),
        s: ci(pc) + (pc.pass ? " — PASS (CI excludes 0)" : " — FAIL") },
      { k: "Program-gene coverage — mouse", v: cov("coverage_mouse"), s: "of the conserved program genes found in each species' matrix" },
      { k: "Program-gene coverage — human", v: cov("coverage_human"), s: "of the conserved program genes found in each species' matrix" },
      { k: "Program-gene coverage — NMR", v: cov("coverage_nmr"), s: "of the conserved program genes found in each species' matrix — scores rest on a partial gene set" }
    ];
    return '<div class="stat-row">' + cards.map(function (c) {
      return '<div class="stat"><div class="k">' + esc(c.k) + '</div><div class="v">' + c.v +
        '</div><div class="s">' + esc(c.s) + "</div></div>";
    }).join("") + "</div>";
  }

  function lifespansHtml() {
    var rows = DATA.signature.lifespans.map(function (r) {
      return "<tr><td>" + esc(r["Unnamed: 0"] || "") + "</td><td>" + esc(r.common_name) +
        '</td><td class="num">' + f(r.max_lifespan_years, 1) + " years</td></tr>";
    }).join("");
    return "<table class='data'><thead><tr><th>Species</th><th>Common name</th><th class='num'>Max lifespan (AnAge)</th></tr></thead><tbody>" +
      rows + "</tbody></table>";
  }

  var SIG_COLORS = {
    "Mouse (whole-body)": "#1f77b4",
    "Mouse": "#9ecae1",
    "Mouse (epidermal)": "#c6dbef",
    "Human": "#d62728",
    "Human (fibroblast)": "#f4a3a3",
    "Naked mole-rat": "#2ca02c"
  };

  function drawSignatureScatter(canvas) {
    var scores = DATA.signature.scores;
    var W = canvas.width, H = canvas.height, padL = 52, padR = 16, padT = 14, padB = 40;
    var ctx = canvas.getContext("2d");
    var xs = scores.map(function (r) { return r.age_normalized; });
    var ys = scores.map(function (r) { return r.program_score; });
    var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
    var yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
    var yPad = (yMax - yMin) * 0.08 || 0.1;
    yMin -= yPad; yMax += yPad;
    function X(v) { return padL + (v - xMin) / (xMax - xMin) * (W - padL - padR); }
    function Y(v) { return H - padB - (v - yMin) / (yMax - yMin) * (H - padT - padB); }

    ctx.clearRect(0, 0, W, H);
    ctx.strokeStyle = "#d9dfe6"; ctx.fillStyle = "#46586a"; ctx.lineWidth = 1;
    ctx.font = "11px sans-serif";
    [0, 0.2, 0.4, 0.6].forEach(function (t) {
      ctx.beginPath(); ctx.moveTo(X(t), padT); ctx.lineTo(X(t), H - padB); ctx.stroke();
      ctx.textAlign = "center"; ctx.fillText(t.toFixed(1), X(t), H - padB + 14);
    });
    for (var i = 0; i <= 4; i++) {
      var v = yMin + (yMax - yMin) * i / 4;
      ctx.beginPath(); ctx.moveTo(padL, Y(v)); ctx.lineTo(W - padR, Y(v)); ctx.stroke();
      ctx.textAlign = "right"; ctx.fillText(v.toFixed(2), padL - 5, Y(v) + 3);
    }
    ctx.textAlign = "center"; ctx.fillText("life-stage-normalized age (fraction of species max lifespan)", (W + padL) / 2, H - 6);

    // OLS reference lines for the three primary cohorts, using the CSV slopes
    // (intercepts derived as a = mean(y) - slope * mean(x) — rendering only, no new numbers shown).
    var stats = byKey(DATA.signature.statistics, "comparison");
    [["Mouse (whole-body)", "Mouse (whole-body) aging slope"],
     ["Human", "Human aging slope"],
     ["Naked mole-rat", "Naked mole-rat aging slope"]].forEach(function (pair) {
      var cohort = scores.filter(function (r) { return r.species === pair[0]; });
      var srow = stats[pair[1]];
      if (!cohort.length || !srow) return;
      var mx = 0, my = 0;
      cohort.forEach(function (r) { mx += r.age_normalized; my += r.program_score; });
      mx /= cohort.length; my /= cohort.length;
      var a = my - srow.slope * mx;
      ctx.strokeStyle = SIG_COLORS[pair[0]]; ctx.lineWidth = 1.8;
      ctx.beginPath();
      ctx.moveTo(X(0), Y(a));
      ctx.lineTo(X(xMax), Y(a + srow.slope * xMax));
      ctx.stroke();
    });

    scores.forEach(function (r) {
      ctx.fillStyle = SIG_COLORS[r.species] || "#666";
      ctx.globalAlpha = 0.85;
      ctx.beginPath();
      ctx.arc(X(r.age_normalized), Y(r.program_score), 3.4, 0, Math.PI * 2);
      ctx.fill();
    });
    ctx.globalAlpha = 1;
  }

  function signatureLegendHtml() {
    var seen = {};
    DATA.signature.scores.forEach(function (r) { seen[r.species] = true; });
    return '<div class="legend">' + Object.keys(seen).map(function (sp) {
      return '<span class="key"><span class="swatch" style="background:' + (SIG_COLORS[sp] || "#666") +
        '"></span>' + esc(sp) + "</span>";
    }).join("") + "</div>";
  }

  function renderSignature() {
    var el = document.getElementById("view-signature");
    var sav = DATA.savina;
    var stats = byKey(DATA.signature.statistics, "comparison");
    var g = byKey(DATA.signature.gates, "gate");
    var mw = stats["Mouse (whole-body) aging slope"];
    var nmr = stats["Naked mole-rat aging slope"];
    var tost = stats["NMR TOST equivalence (delta=0.059 = 25% of mouse slope)"];
    var diff = stats["Mouse (whole-body) vs Naked mole-rat slope difference"];
    var perm = stats["NMR aged-vs-young score difference (exact permutation, sample-level)"];
    var covN = g.coverage_nmr.value;
    var humanDonors = DATA.signature.scores.filter(function (r) { return r.species === "Human"; });
    var curveRows = sav.curve.map(function (r) {
      return "<tr><td class='num'>" + f(r.log2fc_cutoff, 2) + "</td><td class='num'>" + r.n_de_pval05 +
        "</td><td class='num'>" + r.n_de_fdr05 + "</td></tr>";
    }).join("");
    var pcRows = sav.positive_controls.map(function (r) {
      return "<tr><td><em>" + esc(r.gene) + "</em></td><td class='num'>" + f(r.log2FC, 3) +
        "</td><td class='num'>" + pfmt(r.pvalue) + "</td><td class='num'>" + pfmt(r.FDR) +
        "</td><td>" + esc(r.direction) + "</td></tr>";
    }).join("");

    el.innerHTML =
      '<h2>1 &middot; Signature contrast: does the conserved aging program track lifespan across species?</h2>' +
      '<p class="lede">The conserved program (' + DATA.program_genes.counts.n_program_up + " up-regulated + " +
      DATA.program_genes.counts.n_program_down + " down-regulated genes with age) scored per sample against " +
      "life-stage-normalized age. <strong>Primary result:</strong> the mouse whole-body cohort (n&nbsp;=&nbsp;" + mw.n +
      " mice) shows a real aging signal (p&nbsp;=&nbsp;" + pfmt(mw.p) + "). " +
      "The human arm is low-powered (n&nbsp;=&nbsp;" + humanDonors.length + " donors) and the naked-mole-rat program-gene coverage is only " +
      (covN * 100).toFixed(1) + "% — and the pre-registered " +
      "TOST equivalence test did <strong>not</strong> establish NMR&ndash;mouse equivalence.</p>" +
      '<div class="callout warn"><strong>Equivalence not established — read the p-values.</strong> ' +
      "NMR slope " + f(nmr.slope, 3) + " (95% CI " + ci(nmr) + ", p&nbsp;=&nbsp;" + pfmt(nmr.p) +
      "); mouse-vs-NMR slope difference " + f(diff.slope, 3) + " (" + ci(diff) + ", p&nbsp;=&nbsp;" + pfmt(diff.p) +
      "); TOST at &delta;&nbsp;=&nbsp;" + f(g.tost_delta.value, 3) + " (25% of the mouse slope) p&nbsp;=&nbsp;" + pfmt(tost.p) +
      "; exact permutation p&nbsp;=&nbsp;" + pfmt(perm.p) +
      ". Absence of a significant NMR slope is <em>not</em> evidence of flatness.</div>" +
      '<div class="grid-2">' +
      '<div class="card"><h3>Signature figure (M6)</h3><figure>' +
      '<img src="figures/m6-signature-figure.png" alt="Conserved program score vs normalized age in mouse, human, and naked mole-rat">' +
      "<figcaption>From <code>figures/m6-signature-figure.png</code> (M6 pipeline output).</figcaption>" +
      "</figure>" + prov("figures/m6-signature-figure.png") + "</div>" +
      '<div class="card"><h3>Regression &amp; equivalence statistics (all rows)</h3>' + statsTableHtml() +
      prov("signature") + "</div>" +
      "</div>" +
      "<h3>M6 gates</h3>" + gatesHtml() +
      '<div class="grid-2">' +
      '<div class="card"><h3>AnAge lifespan denominators</h3>' + lifespansHtml() + prov("signature") + "</div>" +
      '<div class="card"><h3>All ' + DATA.signature.scores.length + ' sample-level scores (CSV slopes overlaid)</h3>' +
      '<canvas id="sig-scatter" width="640" height="420"></canvas>' + signatureLegendHtml() +
      '<p class="provenance">Points: <code>results/m6-lifespan-program-scores.csv</code>; line slopes: <code>results/m6-statistical-summary.csv</code>.</p>' +
      "</div>" +
      "</div>" +
      '<div class="card"><h3>Savina 2022 reproduction — NMR epidermis DE vs |log2FC| cutoff</h3>' +
      '<div class="grid-2"><figure><img src="figures/m6-savina-reproduction.png" alt="NMR DE-count vs log2FC cutoff curve">' +
      "<figcaption>From <code>figures/m6-savina-reproduction.png</code>.</figcaption></figure>" +
      '<div><p>The paper&rsquo;s headline — essentially no global aged-vs-young DE in naked-mole-rat epidermis — is directionally ' +
      "reproduced at stringent cutoffs, but the curve shows DE counts are cutoff-dependent (" +
      sav.curve[0].n_de_pval05 + " genes at |log2FC| &ge; " + f(sav.curve[0].log2fc_cutoff, 2) +
      " at p&lt;0.05, vs " + sav.curve[sav.curve.length - 1].n_de_pval05 + " at |log2FC| &ge; " +
      f(sav.curve[sav.curve.length - 1].log2fc_cutoff, 2) +
      "), so &ldquo;age-resistant&rdquo; is a statement about effect sizes, not a binary fact.</p>" +
      "<figure><img src=\"figures/m6-positive-control-skin.png\" alt=\"Mouse positive control\" style=\"max-width:340px\">" +
      "<figcaption>Positive control: mouse whole-body program score vs age (gate result in the M6 gates panel above).</figcaption></figure></div></div>" +
      "<h4>Canonical aging-control genes in the NMR contrast</h4>" +
      "<table class='data'><thead><tr><th>Gene</th><th class='num'>log2FC</th><th class='num'>p</th><th class='num'>FDR</th><th>Direction</th></tr></thead><tbody>" +
      pcRows + "</tbody></table>" +
      '<div class="callout info" style="margin-top:10px">' + esc(sav.note) + "</div>" +
      prov("savina") + "</div>";

    drawSignatureScatter(document.getElementById("sig-scatter"));
  }

  /* ================================================================
   * VIEW 2 — GENE BROWSER
   * ================================================================ */
  function geneIndex() {
    var set = {};
    DATA.gene_programs.unit_gene_rows.forEach(function (r) { set[r.gene.toLowerCase()] = r.gene; });
    DATA.program_genes.program_up.forEach(function (g) { set[g.toLowerCase()] = g; });
    DATA.program_genes.program_down.forEach(function (g) { set[g.toLowerCase()] = g; });
    DATA.program_genes.conserved_up.forEach(function (g) { set[g.toLowerCase()] = g; });
    DATA.program_genes.conserved_down.forEach(function (g) { set[g.toLowerCase()] = g; });
    return set;
  }

  function membershipBadges(geneLower) {
    var pg = DATA.program_genes;
    var b = [];
    function has(list) { return list.some(function (g) { return g.toLowerCase() === geneLower; }); }
    if (has(pg.conserved_up)) b.push('<span class="chip ok">conserved up-program</span>');
    else if (has(pg.program_up)) b.push('<span class="chip neutral">up with age (program)</span>');
    if (has(pg.conserved_down)) b.push('<span class="chip ok">conserved down-program</span>');
    else if (has(pg.program_down)) b.push('<span class="chip neutral">down with age (program)</span>');
    return b.length ? b.join(" ") : '<span class="chip warn">not a program gene — unit-level DE only</span>';
  }

  function renderGenes() {
    var el = document.getElementById("view-genes");
    var idx = geneIndex();
    var genes = Object.keys(idx).sort().map(function (k) { return idx[k]; });
    var catalog = byKey(DATA.gene_programs.catalog, "unit");

    el.innerHTML =
      '<h2>2 &middot; Gene-level browser (M5 program genes across ' + DATA.gene_programs.catalog.length +
      ' powered mouse units)</h2>' +
      '<p class="lede">Search a gene to see every powered unit in which it is differentially expressed with age ' +
      "(log2FC, FDR; log2FC convention: old vs young, positive = up with age), whether it belongs to the conserved " +
      "aging program, and the unit's program size and ribo/mito fraction.</p>" +
      '<div class="callout info">Program composition (from <code>results/m6-program-genes.csv</code>): ' +
      DATA.program_genes.counts.n_program_up + " up-regulated + " + DATA.program_genes.counts.n_program_down +
      " down-regulated genes; strictly conserved subsets: " + DATA.program_genes.counts.n_conserved_up + " up, " +
      DATA.program_genes.counts.n_conserved_down + " down.</div>" +
      '<div class="gene-search"><input id="gene-input" list="gene-list" placeholder="Type a gene symbol, e.g. Ccl2" autocomplete="off">' +
      '<datalist id="gene-list">' + genes.map(function (g) { return "<option value=\"" + esc(g) + "\">"; }).join("") + "</datalist>" +
      '<span class="hint">' + genes.length + ' searchable symbols (program genes + unit-level DE genes)</span></div>' +
      '<div id="gene-result"></div>';

    var input = document.getElementById("gene-input");
    var result = document.getElementById("gene-result");

    function show(gene) {
      var key = gene.toLowerCase();
      if (!idx[key]) {
        result.innerHTML = '<div class="callout warn">No data for <strong>' + esc(gene) + "</strong> in " +
          "<code>results/m5-gene-level-programs.csv</code> or the program gene lists.</div>";
        return;
      }
      var canonical = idx[key];
      var rows = DATA.gene_programs.unit_gene_rows.filter(function (r) { return r.gene.toLowerCase() === key; });
      var table;
      if (rows.length) {
        table = '<div class="table-wrap"><table class="data"><thead><tr>' +
          "<th>Unit</th><th class='num'>log2FC</th><th class='num'>FDR</th><th class='num'>Unit n_DE</th>" +
          "<th class='num'>Unit frac_ribo_mito</th><th>Ribo/mito gene?</th></tr></thead><tbody>" +
          rows.map(function (r) {
            var c = catalog[r.unit] || {};
            return "<tr><td>" + esc(r.unit.replace("::", " &middot; ").replace(/_/g, " ")) + "</td>" +
              '<td class="num">' + f(r.log2FC, 3) + "</td>" +
              '<td class="num">' + pfmt(r.FDR) + "</td>" +
              '<td class="num">' + (c.n_de == null ? "—" : c.n_de) + "</td>" +
              '<td class="num">' + f(c.frac_ribo_mito, 3) + "</td>" +
              "<td>" + (r.is_ribo_mito ? '<span class="chip warn">yes</span>' : "no") + "</td></tr>";
          }).join("") + "</tbody></table></div>";
      } else {
        table = "<p>This gene is in the program gene lists but has no unit-level DE row in " +
          "<code>results/m5-gene-level-programs.csv</code>.</p>";
      }
      result.innerHTML =
        '<div class="card"><h3><em>' + esc(canonical) + "</em> &nbsp;" + membershipBadges(key) + "</h3>" +
        table + prov("gene_programs") + prov("program_genes") + "</div>";
    }

    input.addEventListener("change", function () { show(input.value.trim()); });
    input.addEventListener("input", function () {
      var v = input.value.trim();
      if (idx[v.toLowerCase()]) show(v);
    });
    show("Ccl2"); // a well-known program gene present in the CSV
  }

  /* ================================================================
   * VIEW 3 — AGING-SIGNAL RANKING
   * ================================================================ */
  function rankingChartSvg() {
    var rows = DATA.ranking.ranked;
    var rowH = 26, top = 8, bottom = 30, left = 250, right = 76;
    var W = left + 560 + right, H = top + rows.length * rowH + bottom;
    var xMax = 1500;
    function X(v) { return left + v / xMax * 560; }
    var tissueColor = { Lung: "#1f77b4", Kidney: "#9467bd" };
    var s = "";
    [0, 250, 500, 1000, 1500].forEach(function (t) {
      s += '<line x1="' + X(t) + '" y1="' + top + '" x2="' + X(t) + '" y2="' + (H - bottom) + '" stroke="#d9dfe6"/>' +
        '<text x="' + X(t) + '" y="' + (H - bottom + 14) + '" text-anchor="middle">' + t + "</text>";
    });
    s += '<text x="' + (left + 280) + '" y="' + (H - 4) + '" text-anchor="middle">pseudobulk DE genes (FDR &lt; 0.05, |log2FC| ≥ 0.25), 95% bootstrap CI</text>';
    rows.forEach(function (r, i) {
      var y = top + i * rowH + rowH / 2;
      var c = tissueColor[r.tissue] || "#666";
      // CI whisker
      s += '<line x1="' + X(Math.min(r.n_de_ci_low, xMax)) + '" y1="' + y + '" x2="' + X(Math.min(r.n_de_ci_high, xMax)) +
        '" y2="' + y + '" stroke="' + c + '" stroke-width="1.4"/>';
      s += '<line x1="' + X(Math.min(r.n_de_ci_low, xMax)) + '" y1="' + (y - 4) + '" x2="' + X(Math.min(r.n_de_ci_low, xMax)) +
        '" y2="' + (y + 4) + '" stroke="' + c + '"/>';
      s += '<line x1="' + X(Math.min(r.n_de_ci_high, xMax)) + '" y1="' + (y - 4) + '" x2="' + X(Math.min(r.n_de_ci_high, xMax)) +
        '" y2="' + (y + 4) + '" stroke="' + c + '"/>';
      // bar
      s += '<rect x="' + left + '" y="' + (y - 8) + '" width="' + Math.max(1, X(r.n_de_genes) - left) +
        '" height="16" fill="' + c + '" fill-opacity="0.85"/>';
      // labels
      s += '<text x="' + (left - 6) + '" y="' + (y + 3) + '" text-anchor="end">' +
        esc(r.cell_ontology_class) + " (" + esc(r.tissue) + ", n=" + r.n_mice_total + " mice)</text>";
      s += '<text x="' + (X(r.n_de_genes) + 5) + '" y="' + (y + 3) + '">' + r.n_de_genes + "</text>";
    });
    return '<div class="chart-holder"><svg width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + " " + H + '">' + s +
      "</svg></div>" +
      '<div class="legend"><span class="key"><span class="swatch" style="background:#1f77b4"></span>Lung</span>' +
      '<span class="key"><span class="swatch" style="background:#9467bd"></span>Kidney</span>' +
      '<span class="key">whiskers = 95% bootstrap CI over mice — wide intervals reflect small n; read ranking, not exact counts</span></div>';
  }

  function renderRanking() {
    var el = document.getElementById("view-ranking");
    var trio = DATA.ranking.calibration_trio[0];
    var sens = DATA.ranking.sensitivity[0];
    var pcRows = DATA.savina.positive_controls.map(function (r) {
      return "<tr><td><em>" + esc(r.gene) + "</em></td><td class='num'>" + f(r.log2FC, 3) + "</td><td class='num'>" +
        pfmt(r.pvalue) + "</td><td class='num'>" + pfmt(r.FDR) + "</td><td>" + esc(r.direction) + "</td></tr>";
    }).join("");
    var tableRows = DATA.ranking.ranked.map(function (r) {
      return "<tr><td>" + esc(r.unit) + "</td><td class='num'>" + r.n_mice_total + "</td><td class='num'>" + r.n_cells_total +
        "</td><td class='num'>" + r.n_de_genes + "</td><td class='num'>" + f(r.n_de_ci_low, 1) + "</td><td class='num'>" +
        f(r.n_de_ci_high, 1) + "</td><td class='num'>" + f(r.median_abs_log2fc, 3) + "</td></tr>";
    }).join("");

    el.innerHTML =
      "<h2>3 &middot; Mouse aging-signal ranking (M2)</h2>" +
      '<p class="lede">' + DATA.ranking.n_units_powered + " of " + DATA.ranking.n_units_evaluated +
      " tissue &times; cell-type units passed the eligibility rule (" + esc(DATA.ranking.eligibility_rule) +
      "). Ranking metric: pseudobulk DE gene count on the 3m vs 18m contrast, with bootstrap-over-mice CIs.</p>" +
      '<div class="stat-row">' +
      '<div class="stat"><div class="k">Units evaluated</div><div class="v">' + DATA.ranking.n_units_evaluated + '</div><div class="s">all tissue × cell-type combos in the droplet arm</div></div>' +
      '<div class="stat"><div class="k">Units powered</div><div class="v">' + DATA.ranking.n_units_powered + '</div><div class="s">≥3 mice and ≥50 cells per group</div></div>' +
      '<div class="stat"><div class="k">Calibration trio — top unit</div><div class="v">' + esc(trio.top_unit) + '</div><div class="s">observed ' + trio.observed_aging_de +
        " DE vs label-permutation null " + f(trio.permutation_null_mean_de, 1) + " vs sex-effect floor (3m) " + trio.sex_floor_3m_de + " DE</div></div>" +
      "</div>" +
      '<div class="card"><h3>Ranking of the ' + DATA.ranking.n_units_powered + ' powered units</h3>' + rankingChartSvg() + "</div>" +
      '<div class="card"><h3>Full ranking table</h3><div class="table-wrap"><table class="data"><thead><tr>' +
      "<th>Unit</th><th class='num'>Mice</th><th class='num'>Cells</th><th class='num'>n_DE</th>" +
      "<th class='num'>CI low</th><th class='num'>CI high</th><th class='num'>median |log2FC|</th>" +
      "</tr></thead><tbody>" + tableRows + "</tbody></table></div>" + prov("ranking") + "</div>" +
      '<div class="grid-2">' +
      '<div class="card"><h3>Sensitivity runs</h3><table class="data"><thead><tr><th>Run</th><th class="num">DE genes</th></tr></thead><tbody>' +
      "<tr><td>Primary 3m vs 18m (" + esc(sens.unit) + ")</td><td class='num'>" + sens.primary_3m_vs_18m_de + "</td></tr>" +
      "<tr><td>Extended bins (1m+3m vs 18m+21m+24m+30m)</td><td class='num'>" + sens.extended_bins_de + "</td></tr>" +
      "<tr><td>Male-only 18m vs 24m (" + esc(sens.male_only_unit) + ")</td><td class='num'>" + sens.male_only_18m_vs_24m_de + "</td></tr>" +
      "</tbody></table><p class='provenance'>Source: <code>results/m2-sensitivity-summary.csv</code>.</p></div>" +
      '<div class="card"><h3>Positive controls (canonical aging genes)</h3>' +
      "<table class='data'><thead><tr><th>Gene</th><th class='num'>log2FC</th><th class='num'>p</th><th class='num'>FDR</th><th>Direction</th></tr></thead><tbody>" +
      pcRows + "</tbody></table>" +
      '<p class="provenance">Source: <code>results/m2-calibration-positive-controls.csv</code>. Note: Igfbp3 does not replicate (FDR ' +
      f(DATA.savina.positive_controls[0].FDR, 3) + ") — see Caveats.</p></div>" +
      "</div>";
  }

  /* ================================================================
   * VIEW 4 — CONSERVATION / PRIOR WORK
   * ================================================================ */
  function renderConservation() {
    var el = document.getElementById("view-conservation");
    var vc = DATA.agreement.verdict_counts;
    var agreeRows = DATA.agreement.agreement.map(function (r) {
      var chip = r.verdict === "agree" ? '<span class="chip ok">agree</span>'
        : r.verdict === "partial" ? '<span class="chip partial">partial</span>'
        : '<span class="chip bad">disagree</span>';
      return "<tr><td>" + esc(r.unit) + "</td><td>" + esc(r.paper_column) + "</td><td class='num'>" + r.n_shared +
        "</td><td class='num'>" + r.n_paper_sig_bh0.01 + "</td><td class='num'>" + f(r.spearman_rho, 3) +
        "</td><td class='num'>" + f(r.sign_concordance_paper_sig, 3) + "</td><td class='num'>" +
        f(r.our_nominal_recovery_in_paper_sig, 3) + "</td><td>" + chip + "</td><td>" + esc(r.reason) + "</td></tr>";
    }).join("");
    var divRows = DATA.agreement.divergence.map(function (r) {
      return "<tr><td>" + esc(r.unit.replace(/_/g, " ")) + "</td><td class='num'>" + r.n_de + "</td><td class='num'>" + r.n_up +
        "</td><td class='num'>" + r.n_down + "</td><td class='num'>" + f(r.median_abs_log2fc, 3) + "</td><td class='num'>" +
        f(r.frac_ribo_mito, 3) + "</td><td>" + '<span class="chip warn">human arm untestable</span>' + "</td></tr>";
    }).join("");
    var hpRows = DATA.agreement.human_power_grid.map(function (r) {
      return "<tr><td class='num'>" + r.n_donors + "</td><td class='num'>" + r.donor_age_sd_years + "</td><td class='num'>" +
        r.biological_dispersion + "</td><td class='num'>" + f(r.se_slope_per_year, 4) + "</td><td class='num'>" +
        f(r.min_detectable_slope_per_year, 4) + "</td><td class='num'>" + f(r.min_detectable_30y_log2fc, 3) + "</td><td>" +
        esc(r.power_status) + "</td></tr>";
    }).join("");

    el.innerHTML =
      "<h2>4 &middot; Conservation &amp; agreement with prior work (M3 / M5)</h2>" +
      '<p class="lede">Agreement between our DE results and the Tabula Muris Senis paper&rsquo;s published aging coefficients, ' +
      "and the cross-species conservation table. The human arm of the conservation analysis is reported as a " +
      "<strong>first-class result: untestable with obtainable public data</strong> — not an apology.</p>" +
      '<div class="callout bad"><strong>The human arm is untestable — a result, not a gap in effort.</strong> ' +
      "Every unit in the divergence table carries the status &ldquo;" +
      esc(DATA.agreement.divergence[0].status) + "&rdquo;. The Tabula Sapiens cohort is adult with no usable age " +
      "contrast, the SenNet public search returned no qualifying dataset (documented in M4), and the M4 power grid " +
      "below flags every human scenario as low-powered. Therefore no conservation/divergence call is made for any unit.</div>" +
      '<div class="stat-row">' +
      '<div class="stat"><div class="k">Verdicts vs the paper&rsquo;s published coefficients</div><div class="v">' +
        vc.agree + " agree &middot; " + vc.partial + " partial &middot; " + (vc.disagree || 0) + " disagree</div>" +
        '<div class="s">across the ' + DATA.agreement.agreement.length +
        ' powered units; sign concordance on the paper&rsquo;s FDR&lt;0.01 genes</div></div>' +
      "</div>" +
      '<div class="card"><h3>Agreement table (M3)</h3><div class="table-wrap"><table class="data"><thead><tr>' +
      "<th>Unit</th><th>Paper column</th><th class='num'>Shared genes</th><th class='num'>Paper sig (BH 0.01)</th>" +
      "<th class='num'>Spearman &rho;</th><th class='num'>Sign concordance</th><th class='num'>Nominal recovery</th>" +
      "<th>Verdict</th><th>Method note</th></tr></thead><tbody>" + agreeRows + "</tbody></table></div>" +
      prov("agreement") + "</div>" +
      '<div class="card"><h3>Divergence table (M5) — mouse catalogued; human arm untestable for all ' +
      DATA.agreement.divergence.length + ' units</h3>' +
      '<div class="table-wrap"><table class="data"><thead><tr><th>Unit</th><th class="num">n_DE</th><th class="num">n_up</th>' +
      "<th class='num'>n_down</th><th class='num'>median |log2FC|</th><th class='num'>frac_ribo_mito</th><th>Status</th>" +
      "</tr></thead><tbody>" + divRows + "</tbody></table></div>" + prov("agreement") + "</div>" +
      '<div class="grid-2">' +
      '<div class="card"><h3>Human power grid (M4)</h3>' +
      '<div class="table-wrap"><table class="data"><thead><tr><th class="num">n donors</th><th class="num">Age SD (y)</th>' +
      "<th class='num'>Dispersion</th><th class='num'>SE slope/y</th><th class='num'>Min detectable slope/y</th>" +
      "<th class='num'>Min detectable 30y log2FC</th><th>Power status</th></tr></thead><tbody>" + hpRows +
      "</tbody></table></div>" + prov("agreement") + "</div>" +
      '<div class="card"><h3>Program recurrence (M5)</h3><figure>' +
      '<img src="figures/m5-program-recurrence.png" alt="Program gene recurrence across the powered units">' +
      "<figcaption>From <code>figures/m5-program-recurrence.png</code>. " +
      '(<code>figures/m5-housekeeping-guard.png</code> is not present in this repo state, so the housekeeping-guard figure is not shown.)</figcaption>' +
      "</figure>" + prov("figures/m5-program-recurrence.png") + "</div>" +
      "</div>";
  }

  /* ================================================================
   * VIEW 5 — UMAP EXPLORER
   * ================================================================ */
  var UMAP_PALETTES = {
    age: ["#2b6cb0", "#4299e1", "#a0cbe8", "#e8a04c", "#d97b29", "#b3541e"],
    tissue: ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
             "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
             "#98df8a", "#ff9896", "#c5b0d5", "#c49c94"],
    sex: ["#d62728", "#1f77b4"],
    cell_type: null // hashed hues
  };
  function hashHue(i) { return "hsl(" + ((i * 47) % 360) + ",62%,48%)"; }

  function drawUmap(canvas, mode) {
    var um = DATA.umap_subsample;
    var W = canvas.width, H = canvas.height, pad = 14;
    var ctx = canvas.getContext("2d");
    var xs = um.x, ys = um.y;
    var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
    var yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
    function X(v) { return pad + (v - xMin) / (xMax - xMin) * (W - 2 * pad); }
    function Y(v) { return pad + (v - yMin) / (yMax - yMin) * (H - 2 * pad); }
    ctx.clearRect(0, 0, W, H);
    var codes = um[mode === "cell_type" ? "cell_type" : mode];
    var palette = UMAP_PALETTES[mode];
    for (var i = 0; i < xs.length; i++) {
      var c = codes[i];
      ctx.fillStyle = palette ? palette[c % palette.length] : hashHue(c);
      ctx.globalAlpha = 0.72;
      ctx.beginPath();
      ctx.arc(X(xs[i]), Y(ys[i]), mode === "tissue" || mode === "cell_type" ? 1.8 : 2.4, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
  }

  function umapLegendHtml(mode) {
    var um = DATA.umap_subsample;
    var cats = um.categories[mode === "cell_type" ? "cell_type" : mode];
    var palette = UMAP_PALETTES[mode];
    var items = cats.map(function (name, i) {
      var col = palette ? palette[i % palette.length] : hashHue(i);
      return '<span class="key"><span class="swatch" style="background:' + col + '"></span>' + esc(name) + "</span>";
    });
    return '<div class="legend">' + items.join("") + "</div>";
  }

  function renderUmap() {
    var el = document.getElementById("view-umap");
    var um = DATA.umap_subsample;
    var sa = DATA.confound.sex_age;
    var saRows = sa.map(function (r) {
      var singleSex = (r.cells_female === 0 || r.cells_male === 0) && r.age_group !== "All";
      return "<tr" + (singleSex ? ' class="hl"' : "") + "><td>" + esc(r.age_group) + "</td><td class='num'>" + r.cells_female +
        "</td><td class='num'>" + r.cells_male + "</td><td class='num'>" + r.cells_total + "</td><td class='num'>" + r.mice_female +
        "</td><td class='num'>" + r.mice_male + "</td><td class='num'>" + r.mice_total + "</td></tr>";
    }).join("");
    var compRows = DATA.confound.composition.map(function (r) {
      return "<tr><td>" + esc(r.cell_ontology_class) + "</td><td class='num'>" + f(r.fraction_paper_raw, 4) +
        "</td><td class='num'>" + f(r.fraction_our_pipeline, 4) + "</td><td class='num'>" + f(r.absolute_diff, 4) + "</td></tr>";
    }).join("");
    var pt = um.meta.per_tissue_counts;
    var ptText = Object.keys(pt).map(function (t) { return esc(t) + " " + pt[t]; }).join(" · ");
    var singleSex = DATA.confound.sex_age.filter(function (r) {
      return r.age_group !== "All" && (r.cells_female === 0 || r.cells_male === 0);
    }).map(function (r) { return r.age_group; });
    var mixed = DATA.confound.sex_age.filter(function (r) {
      return r.age_group !== "All" && r.cells_female > 0 && r.cells_male > 0;
    }).map(function (r) { return r.age_group; });

    el.innerHTML =
      "<h2>5 &middot; Atlas UMAP explorer (M1)</h2>" +
      '<p class="lede">Seeded per-tissue proportional subsample of the full QC&rsquo;d atlas (' +
      um.meta.n_cells_total_in_object.toLocaleString("en-US") + " cells): " +
      um.meta.n_sampled + " cells sampled (seed " + um.meta.seed +
      ", without replacement). Coordinates are the M1 <code>X_umap</code> rounded to 2 decimals for display.</p>" +
      '<div class="callout warn"><strong>Sex &times; age confounding is structural in the droplet arm.</strong> ' +
      "Age groups " + singleSex.join(", ") + " are single-sex; only " + mixed.join(", ") +
      " is/are sex-balanced (highlighted below). Primary DE conclusions use the 3m vs 18m contrast; extended-bin runs are sensitivity analyses.</div>" +
      '<div class="card"><div class="umap-controls"><label>Color by: ' +
      '<select id="umap-mode"><option value="age" selected>Age group</option><option value="tissue">Tissue</option>' +
      '<option value="sex">Sex</option><option value="cell_type">Cell type</option></select></label>' +
      '<span class="hint" style="color:var(--ink-soft);font-size:12.5px">sampling: ' + esc(um.meta.sampling) + "</span></div>" +
      '<canvas id="umap-canvas" width="820" height="560"></canvas>' +
      '<div id="umap-legend"></div>' +
      '<p class="provenance">Per-tissue sample sizes: ' + ptText + ". Source: <code>data/interim/m1-processed.h5ad</code> " +
      '(obsm[\'X_umap\'] + obs), recorded in <code>data/manifest.json</code> with SHA-256.</p></div>' +
      '<div class="grid-2">' +
      '<div class="card"><h3>Sex &times; age table (cells and mice per group)</h3>' +
      "<table class='data'><thead><tr><th>Age group</th><th class='num'>Cells F</th><th class='num'>Cells M</th>" +
      "<th class='num'>Cells total</th><th class='num'>Mice F</th><th class='num'>Mice M</th><th class='num'>Mice total</th></tr></thead><tbody>" +
      saRows + "</tbody></table>" + prov("confound") + "</div>" +
      '<div class="card"><h3>Cell-type composition vs paper</h3>' +
      '<details><summary>Show ' + DATA.confound.composition.length + " cell types (fractions paper vs this pipeline)</summary>" +
      '<div class="table-wrap" style="margin-top:8px"><table class="data"><thead><tr><th>Cell type</th>' +
      "<th class='num'>Fraction (paper)</th><th class='num'>Fraction (ours)</th><th class='num'>|diff|</th></tr></thead><tbody>" +
      compRows + "</tbody></table></div></details>" + prov("confound") + "</div>" +
      "</div>";

    var canvas = document.getElementById("umap-canvas");
    var select = document.getElementById("umap-mode");
    var legend = document.getElementById("umap-legend");
    function refresh() {
      var mode = select.value;
      drawUmap(canvas, mode);
      legend.innerHTML = umapLegendHtml(mode);
    }
    select.addEventListener("change", refresh);
    refresh();
  }

  /* ================================================================
   * CAVEATS MODAL
   * ================================================================ */
  function caveat(num, text, srcs) {
    return "<li><strong>" + num + ".</strong> " + text +
      '<br><span class="src">Source: ' + srcs.map(function (s) { return "<code>" + esc(s) + "</code>"; }).join(", ") + "</span></li>";
  }

  function renderCaveats() {
    var g = byKey(DATA.signature.gates, "gate");
    var stats = byKey(DATA.signature.statistics, "comparison");
    var sig = DATA.signature;
    var humSlope = stats["Human aging slope"], fibSlope = stats["Human (fibroblast) aging slope"];
    var mSkin = stats["Mouse aging slope"];
    var tost = stats["NMR TOST equivalence (delta=0.059 = 25% of mouse slope)"];
    var perm = stats["NMR aged-vs-young score difference (exact permutation, sample-level)"];
    var diff = stats["Mouse (whole-body) vs Naked mole-rat slope difference"];
    var covN = g.coverage_nmr.value, covH = g.coverage_human.value;
    var humanDonors = sig.scores.filter(function (r) { return r.species === "Human"; });
    var humanAges = humanDonors.map(function (r) { return r.age_years; });
    var nmr = sig.scores.filter(function (r) { return r.species === "Naked mole-rat"; });
    var nmrAges = {};
    nmr.forEach(function (r) { nmrAges[r.age_years] = (nmrAges[r.age_years] || 0) + 1; });
    var nmrAgeVals = Object.keys(nmrAges).map(Number);
    var ageStr = nmrAgeVals.sort(function (a, b) { return a - b; }).map(function (a) { return a + "y ×" + nmrAges[a]; }).join(", ");
    var nmrLifespan = sig.lifespans.filter(function (r) { return r.common_name === "Naked mole-rat"; })[0].max_lifespan_years;
    var nmrAgedFrac = (Math.max.apply(null, nmrAgeVals) / nmrLifespan * 100).toFixed(0);
    var singleSex = DATA.confound.sex_age.filter(function (r) {
      return r.age_group !== "All" && (r.cells_female === 0 || r.cells_male === 0);
    }).map(function (r) { return r.age_group; });

    var items = [
      caveat(1, "<strong>Human arm untestable with obtainable public data.</strong> All 17 divergence-table units carry the status “" +
        DATA.agreement.divergence[0].status + "”; the M4 power grid flags every scenario “" +
        DATA.agreement.human_power_grid[0].power_status + "”. No conservation/divergence call is possible or made.",
        ["results/m5-divergence-table.csv", "results/m4-human-power.csv"]),
      caveat(2, "<strong>Human cohort is tiny.</strong> n = " + humanDonors.length + " TS donors (ages " + Math.min.apply(null, humanAges) +
        "–" + Math.max.apply(null, humanAges) + "y, all adult — no young pole). Human slopes are uninformative: all-cells " +
        f(humSlope.slope, 3) + " (" + ci(humSlope) + ", p = " + pfmt(humSlope.p) + "); fibroblast-only " + f(fibSlope.slope, 3) +
        " (" + ci(fibSlope) + ", p = " + pfmt(fibSlope.p) + ").",
        ["results/m6-lifespan-program-scores.csv", "results/m6-statistical-summary.csv"]),
      caveat(3, "<strong>NMR program coverage is only " + (covN * 100).toFixed(1) + "%</strong> (human " + (covH * 100).toFixed(1) +
        "%, mouse " + (g.coverage_mouse.value * 100).toFixed(1) + "%). NMR program scores rest on a partial gene set and are not directly comparable to mouse scores.",
        ["results/m6-gates.csv"]),
      caveat(4, "<strong>NMR–mouse equivalence is NOT established.</strong> The pre-registered TOST at δ = " + f(g.tost_delta.value, 4) +
        " (25% of the mouse slope) gives p = " + pfmt(tost.p) + "; the slope difference is " + f(diff.slope, 3) + " (" + ci(diff) +
        ", p = " + pfmt(diff.p) + "). Failure to reject a zero NMR slope is absence of evidence, not evidence of flatness.",
        ["results/m6-statistical-summary.csv", "results/m6-gates.csv"]),
      caveat(5, "<strong>No detectable NMR aged-vs-young signal at all in these data:</strong> exact sample-level permutation p = " +
        pfmt(perm.p) + " (n = " + perm.n + ").",
        ["results/m6-statistical-summary.csv"]),
      caveat(6, "<strong>NMR ages are group means with only " + nmrAgeVals.length + " distinct values</strong> (" + ageStr + "); the “aged” group (~" +
        nmrAgedFrac + "% of the " + f(nmrLifespan, 1) + "y AnAge max lifespan) is middle-aged in NMR terms, and the sample is cultured skin fibroblasts only — no multi-tissue NMR atlas exists publicly.",
        ["results/m6-lifespan-program-scores.csv", "results/m6-species-lifespans.csv"]),
      caveat(7, "<strong>Igfbp3 — a signature gene and the Savina paper's key NMR claim — is absent from the processed GEO matrix</strong>, so it is untestable in NMR here; in the mouse calibration it is also not DE (FDR = " +
        f(DATA.savina.positive_controls[0].FDR, 3) + ").",
        ["results/m2-calibration-positive-controls.csv", "scripts/run_m6_pipeline.py"]),
      caveat(8, "<strong>Sex × age confounding in the mouse droplet arm:</strong> age groups " + singleSex.join(", ") +
        " are single-sex. Only 3m vs 18m is sex-balanced; everything else is a sensitivity analysis.",
        ["results/m1-sex-age-table.csv"]),
      caveat(9, "<strong>Most sub-cohorts are underpowered.</strong> Only the mouse whole-body cohort is powered (n = " +
        stats["Mouse (whole-body) aging slope"].n + ", slope " +
        f(stats["Mouse (whole-body) aging slope"].slope, 3) + ", p = " + pfmt(stats["Mouse (whole-body) aging slope"].p) +
        "). Mouse skin n = " + stats["Mouse aging slope"].n + " (p = " + pfmt(mSkin.p) +
        "), human n = " + humanDonors.length + " (see #2), NMR n = " + nmr.length +
        " with " + nmrAgeVals.length + " age points.",
        ["results/m6-statistical-summary.csv"]),
      caveat(10, "<strong>Cross-sectional design with survivorship bias:</strong> 24m/30m mice are biological survivors; datasets are cross-sectional cohorts, not longitudinal tracking.",
        ["results/m1-sex-age-table.csv", "notes/age-bins.md"]),
      caveat(11, "<strong>Correlation ≠ causation:</strong> every transcriptomic shift here is observational. Whether program genes drive or respond to aging is untested.",
        ["report/conserved-aging-programs.md (Limitations)"])
    ];
    document.getElementById("caveats-body").innerHTML =
      "<p>These limitations are part of the results. The app shows the numbers anyway, with their intervals and p-values, " +
      "because that is the honest picture.</p><ol class='caveats'>" + items.join("") + "</ol>";
  }

  function wireModal() {
    var modal = document.getElementById("caveats-modal");
    document.getElementById("caveats-btn").addEventListener("click", function () { modal.hidden = false; });
    document.getElementById("caveats-close").addEventListener("click", function () { modal.hidden = true; });
    modal.addEventListener("click", function (e) { if (e.target === modal) modal.hidden = true; });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") modal.hidden = true; });
  }

  /* ---------------- tabs ---------------- */
  function wireTabs() {
    var tabs = Array.prototype.slice.call(document.querySelectorAll(".tab"));
    function activate(name) {
      tabs.forEach(function (t) { t.classList.toggle("active", t.dataset.view === name); });
      document.querySelectorAll(".view").forEach(function (v) { v.hidden = v.id !== "view-" + name; });
      if (history.replaceState) history.replaceState(null, "", "#" + name);
    }
    tabs.forEach(function (t) { t.addEventListener("click", function () { activate(t.dataset.view); }); });
    var initial = (location.hash || "#signature").slice(1);
    if (!document.getElementById("view-" + initial)) initial = "signature";
    activate(initial);
  }

  /* ---------------- boot ---------------- */
  loadAll().then(function () {
    renderSignature();
    renderGenes();
    renderRanking();
    renderConservation();
    renderUmap();
    renderCaveats();
    wireTabs();
    wireModal();
  }).catch(function (err) {
    document.getElementById("load-error").hidden = false;
    if (window.console) console.error(err);
  });
})();
