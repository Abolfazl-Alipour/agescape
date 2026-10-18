# Agescape Webapp: MVP Specification v2 (Graph-First)

**Product name (working):** Agescape, Conserved Aging Explorer
**Version:** Spec v2.0, 2026-10-18. Supersedes v1.0 in full.
**Basis:** v1.0 as revised by the 10-round hostile review in `agescape-webapp-review-10-rounds.md`. Where v1 and v2 conflict, v2 governs.
**Owner:** Reza Alipour (personal project; Kavox Labs affiliation noted, not the product owner)
**Status:** Ready to build, after the Step 0 plot data audit. Public release only after the usability gate (Section 10) and Reza's staging click-through and explicit go-ahead.
**Source of truth:** The `agescape` repository: committed files in `results/` and `figures/`, plus the plot data exports defined in Section 7. No number, point, or distribution may appear in the app that does not trace to one of those files.

---

## 1. The one job

The app has exactly one job: **convert a 60-second visit from a senior scientist at an aging company into a desire to talk to Reza.**

Every design decision is subordinate to that job. Researcher utility (downloads, shareable views, reusable data) is a welcome by-product and is served through secondary controls. It never drives layout. If a feature does not help a smart, busy, sceptical visitor understand and trust the work faster, it is cut or demoted.

The visitor's path is assumed to be: a link from a resume, LinkedIn, or a message, opened on a phone, by someone who owes the page nothing and will leave the moment it asks for effort.

## 2. Design law

These rules are binding. A build that violates one fails acceptance even if every number is right.

1. **Graphs carry the story. Text explains only what a graph cannot.** Each chapter is one claim, one main graph, and supporting graphs. Prose is rationed by the text budget in Section 5.
2. **Graph titles are conclusions, not topics.** Write "Aging signal is powered in only 2 of 16 tissues", never "Eligibility overview". A visitor who reads only graph titles gets the whole story correctly.
3. **Honesty lives in the geometry.** Absence, weakness, and uncertainty are drawn, not footnoted: hollow markers for untestable units, dashed lines and wide bands for weak evidence, a coverage meter inside the chart, an interval that visibly fails to enter the equivalence zone. Words confirm what the eye already saw. Detail and caveats live one tap down in expanders.
4. **One page, one direction.** The story is a single scroll through five chapters in the fixed order of Section 4. There is no navigation decision before the visitor has received value. The only other surface is an overlay for Methods and a slide-over panel for unit detail.
5. **Three gestures only.** Scroll, tap a graph element, tap an expander. No hover-dependent information, no drag interactions, no multi-step flows, no hamburger menu.
6. **Direct labels beat legends.** Label series, values, and landmarks on the graph itself wherever space allows. A legend is a failure of labelling, permitted only when direct labels genuinely cannot fit.
7. **Mobile is the design surface.** Every graph is designed at 360 px wide first and must be fully legible there. Desktop layouts are adaptations of the phone layout, never the reverse.
8. **Same data, same colours, everywhere.** Species, direction, tissue, and eligibility encodings are defined once in Section 6 and never vary between chapters.
9. **Every graph can defend itself.** Each graph has a provenance drawer in the same position: source file, rows or filter used, data freeze date, checksum short form, download of exactly the data drawn, and, where a published figure exists, a toggle to compare against it.
10. **Animation teaches or it is cut.** Permitted: a one-time draw-in on first scroll into view, and transitions when a tap changes a graph's state. Both must respect reduced-motion settings. No decorative motion.

## 3. What changed from v1, in one view

- Six routes become one scrolling page, one Methods overlay, one unit slide-over panel.
- Tables stop being the primary UI anywhere. A table appears only inside the unit panel (gene lists, where a table is genuinely the right form) and inside provenance drawers as the accessible data twin of each graph.
- The hero stops being the signature figure. The atlas eligibility map leads, because it is the densest honest visual in the project and because trust must be earned before the ambiguous chapter.
- The power audit stops being sliders and a number. It becomes a tappable tile map of study designs.
- Honesty annotation (banners, mandated warning lines on 17 pages) is replaced by honesty encoding (Section 2, rule 3), backed by exact wording in expanders.
- A Plot Data Contract is added. Graph-first design needs distributions and per-gene data, not only summary tables, and Section 7 names exactly which exports are required, which are conditional, and what each graph falls back to.
- A usability gate with named tasks, time boxes, and a pass bar is added. Failing it blocks launch.

## 4. The page: five chapters in fixed order

Story logic: where the signal is, why it is not noise, what the program is, whether it travels across species, and what it would take to answer the human question. Trust is earned in chapters 1 and 2 and spent in chapter 4.

### Hero and Chapter 1: Where aging can even be measured

**Claim (chapter title):** "Of 163 cell types tested, aging is powered in 17, all in lung and kidney."

**Main graph, G0 Atlas Map.** One row per tissue (16 rows), one tile per tissue by cell-type unit in that row (163 tiles in total). Powered units are solid filled tiles in the tissue colour; ineligible units are hollow outlines. Tile size is uniform (size must not imply importance). A count at each row end reads, for example, "Lung 9 of 21" powered. The visual fact a visitor must absorb in under 5 seconds: only two rows contain any solid tiles.

**Supporting graph, G1 Ranking Bars.** Horizontal bars for the 17 powered units, sorted by DE gene count, bar length equals DE genes, whiskers show the bootstrap interval where the source provides one, coloured by tissue, direct-labelled with the count at the bar end (55, 53, 36, 35, and so on). Tapping any bar opens the unit panel (Section 8).

**Expander:** "Why so few?" Contains G4, the sex by age design grid (below), in miniature, and two sentences: the eligibility floor was pre-registered (at least 3 mice, at least 50 cells per group), and the restriction is a property of the atlas design, not an analysis choice.

**Data:** `results/m2-unit-eligibility.csv`, `results/m2-aging-signal-ranked.csv`.

### Chapter 2: Is it real, or is it noise?

**Claim:** "The top signal is 12 times its own permutation null, and the pipeline agrees with the original paper."

**Main graph, G2 Calibration.** The observed DE count for the top unit (55) drawn against its own controls on one shared axis. Preferred form, if the export exists (Section 7, conditional export C1): a histogram or strip of the permutation distribution, with the observed value as a tall marker line far right, and the sex-effect floor at 0. Fallback form, from summary data alone: a dot plot on a single number line with three direct-labelled points, sex floor 0, permutation mean 4.6, observed 55. Either form must make the distance between signal and null visible at a glance. One honest line sits in the expander, not the caption: the canonical marker genes in this unit were directionally plausible but not individually significant, and the calibration claim rests on the permutation and sex controls.

**Main graph, G3 Agreement.** Two layers. Unconditional layer: a lollipop chart of the 17 powered units, one dot per unit on the Spearman rho axis (0 to 1), coloured by verdict (agree, partial), with the count summary drawn as a stacked strip above it: 13 agree, 4 partial, 0 disagree. Conditional layer, if export C3 exists: a gene-level scatter for the kidney proximal tubule unit, our log2FC against the TMS paper's age coefficient across 13,019 shared genes, rho 0.866 labelled on the plot. This scatter, if available, becomes the chapter's lead visual, because agreement across 13,000 genes is the most persuasive single image in the project.

**Supporting visual, G4 Design Grid.** A small grid, age groups as columns (1, 3, 18, 21, 24, 30 months), sex as rows, each cell showing the mouse count, with male-only and female-only columns visibly hatched. The 3 and 18 month columns are outlined as "the only fair fight". This graph explains the whole project's central constraint in one look and replaces several paragraphs of v1 text.

**Supporting visual, G5 Pseudotime Guard (secondary).** The existing diffusion pseudotime vs age figure, shown small in the expander, with its one-line meaning: differentiation structure does not confound the age axis (rho 0.016), and no trajectory aging claim is made.

**Data:** `results/m2-calibration-trio-summary.csv`, `results/m2-calibration-positive-controls.csv`, `results/m3-agreement-table.csv`, `results/m1-sex-age-table.csv`, `figures/m3-dpt-vs-age.png`.

### Chapter 3: The program

**Claim:** "261 genes move with age across the powered units: 167 up, 94 down."

**Main graph, G6 Program Bars.** One horizontal stacked bar per powered unit (17 bars), length split into up genes (direction colour A) and down genes (direction colour B), direct-labelled with both counts. Sorted by total. Units with zero DE genes appear as empty tracks with a 0 label; they are part of the result.

**Main graph, G7 Gene by Unit Matrix.** A dot matrix: rows are the top program genes (cap of 40 rows on mobile, ranked by number of units they appear in, full set in the expander download), columns are the 17 units, a filled dot means the gene is in that unit's program, dot colour encodes direction. Purpose: show at a glance whether the program is shared across units or is mostly private to one unit. The honest answer, whichever it is, is the visual.

**Conditional graph, G8 Volcano.** For the top unit, Lung adventitial cell, a volcano plot (log2FC against minus log10 adjusted p, DE genes coloured, thresholds drawn as lines) if and only if per-gene statistics are exported (Section 7, conditional export C2). If the export does not exist, this graph is cut. No substitute visualisation of per-gene significance may be invented.

**Expander:** the housekeeping guard as two numbers with a small paired bar (ribosomal and mitochondrial fraction among DE genes vs the 0.017 background), per the M5 catalog, with no interpretive claim attached.

**Data:** `results/m5-mouse-program-catalog.csv`, `results/m6-program-genes.csv`.

### Chapter 4: Does it travel across species?

**Claim:** "In mouse, the program rises with age. In the naked mole-rat, over the span tested, it does not. Whether that difference is real is not yet settled."

**Main graph, G9 Signature Small Multiples.** Three stacked panels (side by side only on wide screens), one per species, sharing one y axis scale so slopes are honestly comparable by eye. Each panel: sample points, fitted line, 95 percent CI band, slope and n direct-labelled on the panel. Mouse: solid line, slope 0.235, p equals 0.008, n equals 23. Human: dashed line, visibly wide band, only 4 points, panel watermarked LOW-POWERED in the plot area itself. NMR: slope 0.017, n equals 5. Inside each panel, a coverage meter: a thin bar filled to the program-gene coverage (mouse 1.00, human 0.73, NMR 0.31) and labelled as coverage. The NMR panel's x axis annotation states the tested span in years (1.1 to 11.3) beside the normalised axis.

**Main graph, G10 The Honest Interval.** A forest-style interval plot that answers "what can we claim?" geometrically. Three rows: mouse slope with its CI (excludes 0), NMR slope with its CI (spans 0 widely), and the slope difference with its CI (minus 0.013 to 0.449, just touching 0), zero drawn as a bold vertical line. A shaded band around the NMR row shows the equivalence margin (plus and minus delta, delta equals 0.0587): the NMR interval visibly extends beyond the band, and a direct label on the plot reads "equivalence not established". A visitor must be able to read the project's exact epistemic position off this one graph: mouse effect real, NMR effect undetected, difference suggestive, equivalence unproven.

**Supporting graph, G11 Lifespan Context.** Three horizontal bars on a log-free linear year scale: mouse maximum lifespan 4 years, NMR 31, human 122.5 (AnAge). The tested span is shaded within the NMR and mouse bars. Purpose: make "the NMR aged group is middle-aged for its species" visible instead of arguable.

**Supporting graph, G12 Savina Curve.** Line graph, DE gene count (y, log scale permitted if labelled) against absolute log2FC cutoff (x), with two marked points: the paper's stated rule (29 genes where the paper reported 0) and the unit-converted cutoff 1.45 (9 genes, p rule). Both standing caveats (Igfbp3 absent from the source file; QC count 10,251 vs 10,232) sit in this graph's expander.

**Data:** `results/m6-lifespan-program-scores.csv`, `results/m6-statistical-summary.csv`, `results/m6-program-genes.csv`, `results/m6-savina-reproduction.csv`, `results/m6-species-lifespans.csv`.

### Chapter 5: Why the human answer is missing, and what would answer it

**Claim:** "With 4 donors, only enormous effects could be seen. Here is the design that could see real ones."

**Main graph, G13 Power Map.** Three side-by-side tile maps (small multiples by dispersion: low, mid, high values as present in the grid). Each map: x axis donor count, y axis donor-age SD in years, one tile per precomputed grid point, tile colour equals the minimum detectable effect expressed as a 30-year log2FC (single sequential colour scale, dark equals only huge effects detectable, light equals sensitive). Tiles that meet the eligibility floor (at least 6 donors and SD at least 12) are outlined. One tile is marked with a cross and direct label: "Tabula Sapiens skin: 4 donors". Tapping any tile updates a readout bar fixed at the chapter's base: donors, SD, dispersion, minimum detectable 30-year log2FC (the exact grid value), and a verdict chip, LOW-POWERED or MEETS FLOOR, driven by the source grid flag and the floor rule together, with any disagreement resolved to LOW-POWERED with a one-line reason in the expander.

**Rules for G13:** Tiles are discrete. No smoothing, no contour interpolation, no gradient that implies values between grid points. The readout may only ever display a value that exists in `results/m4-human-power.csv`. The two anchor tiles must reproduce the report: 4 donors, SD 8, lowest dispersion gives 1.474; 4 donors, SD 20, lowest dispersion gives 0.59, and both read LOW-POWERED.

**Closing strip.** Under the Power Map, the page ends with the constructive turn, in graph form where possible: the minimum eligible design tile pulses or is outlined on the map, and beside it three short lines: what data would settle the human question, a link to the full report, and the author block (Section 9). The last thing on the page is an invitation to talk to Reza, with contact links. That is the conversion the whole page exists to produce.

**Data:** `results/m4-human-power.csv`.

### Methods overlay

Opened from a persistent "Methods" control and from each graph's provenance drawer. Contents: data sources table (dataset, source, IDs, version, download date, from `data/manifest/`), the SenNet paragraph (6,820 datasets queried 2026-10-14, no donor ages in the public index), pipeline summary in eight lines, the deviations log summary (full mouse program scored, whole-body positive control substituted, both flagged), the three load-bearing limitations with a link to the full 18 in the report, AI assistance disclosure matching the report, the 2026-10-14 audit mentioned in one sentence, citation block with copy button, and licence lines. This overlay is the one place long text is allowed, because a visitor who opens it has asked for it.

## 5. Text budget

Measured on the rendered page, excluding the Methods overlay and expanders:

- Chapter title (the claim): at most 15 words each.
- Supporting sentence under each main graph: at most 25 words, at most one per graph.
- Total visible body text on first load, whole page: at most 700 words.
- Graphs occupy at least 70 percent of the visible area of every chapter at both phone and desktop sizes.
- Any statistic that matters must appear on or inside a graph, not only in text. Any text that repeats a graph's title claim verbatim is deleted.

## 6. Visual system

- **Species colours (fixed):** mouse, a strong blue; human, a muted purple; naked mole-rat, a warm orange. Used in G9, G10, G11, and nowhere else, so species colour always means species.
- **Direction colours (fixed, colour-blind safe pair):** up with age, vermillion; down with age, blue-grey. Used in G6, G7, G8 only. Never mixed with species colours in the same graph.
- **Eligibility encoding (fixed):** solid fill means powered or observed; hollow outline means ineligible, untestable, or not detected. This single convention is taught by G0 and reused in G9 (hollow points are not used; the human panel uses dashed line plus watermark instead) and G13 (outlined tiles meet the floor; the encoding note appears once, in the G13 expander, to avoid contradiction: in G13 the outline means eligibility, stated on the readout).
- **Tissue colours:** a categorical palette assigned per tissue in G0 and G1 and used consistently if tissues recur. Lung and Kidney must be visually distinct from each other and from the 14 tissues with no powered units, whose rows stay in light neutrals so the two active rows carry the colour weight.
- **Typography:** one sans family, large graph titles (they are headlines), numerals in tabular form. No stock imagery, no decorative graphics, no gradient backgrounds behind data.
- **Charts:** one rendering library only, chosen at Step 0 (ECharts or Observable Plot). All charts are generated from the Section 7 exports at build or load time from static JSON. No chart may be a traced copy of a published PNG presented as data; the published PNG appears only in the compare toggle (Section 8).

## 7. Plot Data Contract

### Step 0 audit, expanded

Before any UI work, the builder produces `app/DATA_CONTRACT.md` in the agescape repo recording, for every file below: exact filename, header columns, row count, one sample row, SHA-256. Every value in this spec is then verified against the files. A file wins any disagreement; the spec is corrected and Reza is told before building continues.

### Required exports (must exist for v2 to ship)

Produced by `scripts/export_app_data.py` from committed results only, written to `app/data/` as JSON, deterministic on re-run:

| Export | Feeds | Source |
|---|---|---|
| `atlas.json` | G0, G1 | m2 eligibility plus ranked tables |
| `calibration.json` | G2 fallback form | m2 calibration trio summary |
| `agreement.json` | G3 lollipop layer | m3 agreement table |
| `design-grid.json` | G4 | m1 sex by age table |
| `program.json` | G6, G7 | m5 catalog plus m6 program genes |
| `signature.json` | G9, G10, G11 | m6 scores, statistical summary, lifespans |
| `savina.json` | G12 | m6 Savina reproduction |
| `power-grid.json` | G13 | m4 human power grid, allowed values taken from the file |
| `provenance.json` | All drawers | source filenames, checksums, export date |

### Conditional exports (each unlocks exactly one graph upgrade)

- **C1, permutation replicates** (the per-replicate DE counts behind the 4.6 mean). Unlocks the G2 histogram form. Producing C1 must be a pure export of values the pipeline already computed. If the pipeline did not retain them, G2 ships in fallback form and re-running permutations is a separate analysis decision for Reza, out of scope for the app build.
- **C2, per-gene DE statistics for the top unit** (gene, log2FC, adjusted p). Unlocks G8, the volcano. If no committed file holds these, G8 is cut.
- **C3, paired effect sizes for the kidney proximal tubule unit** (our log2FC and the paper's age coefficient per shared gene). Unlocks the G3 scatter as chapter lead. If unavailable, G3 ships as the lollipop layer alone.

### Hard rules

- No export may be produced by new analysis choices made by the app builder. Exports transcribe or reshape committed results. Anything requiring a new statistical decision is flagged to Reza and deferred.
- Total `app/data/` under 5 MB. G7's gene cap and any sampling for display must be a stated rule in the contract (for example top 40 genes by unit count), never silent dropping, and the full data stays one tap away in the drawer download.
- `scripts/check_app_data.py` asserts every anchor value in Section 10 against the exports. The build fails if the check fails.

## 8. Interaction specification

- **Scroll story.** Chapters appear in the fixed order of Section 4. A slim progress rail (five dots, thumb-reachable on mobile, side rail on desktop) shows position; tapping a dot scrolls to that chapter. This rail and two persistent links (Report, Code) are the entire navigation.
- **Unit panel.** Tapping a G1 bar or a powered G0 tile opens a slide-over panel for that unit: its G6 row enlarged, its gene table (searchable, downloadable; the one true table in the app), its agreement values, its housekeeping pair, and its conservation line: "Human replication: UNTESTABLE with obtainable public data", linking to Chapter 5. The URL hash updates to the unit slug so the view is shareable and the back button closes the panel. Panel content is generated from the same exports; there are not 17 hand-built anythings.
- **Expanders.** Every graph has two consistent controls in the same positions: "Why this holds, and its limits" (expander with the exact statistics, caveats, and approved wording) and a provenance icon opening the drawer. A visitor learns the pattern once, on G0, and it never changes.
- **Provenance drawer.** Source file, rows or filter, freeze date, checksum short form, download of the drawn data as CSV, and, for G2, G9, and G12, a "Published figure" toggle that swaps the interactive graph for the original figure PNG and back. Any visual discrepancy between the two is a build failure to be reported, not a cosmetic detail.
- **Power Map readout.** Tapping a tile is the only stateful interaction on the page besides the panel and expanders. The readout is a fixed bar, not a tooltip, so it works on touch and persists for screenshots.
- **No dead ends.** Every chapter ends by visually handing off to the next (the final element of each chapter is the first colour or motif of the next). The page bottom contains the author block and contact links.

## 9. Author and conversion block

Page footer, in this order: one line stating the project in plain words; Reza's name, role line (Data Scientist at Meta, PhD Neuroscience, PharmD), and links to LinkedIn, Google Scholar, GitHub, the full report, and a contact email; one line noting the Kavox Labs affiliation and that this is a personal research project; data credits in one line; the line "Research tool, not medical advice"; app version and data freeze date (2026-10-14).

No newsletter form, no analytics-driven popups, no calls to action other than reading the report, viewing the code, or contacting Reza.

## 10. Acceptance tests

**Numeric fidelity (automated, build-blocking):** 163 units and 17 powered, powered confined to Lung and Kidney; ranking order and counts 55, 53, 36, 35 at the head; calibration 55, 4.6, 0; agreement 13, 4, 0 and kidney proximal tubule rho 0.866; program 167 up, 94 down, 312 gene-unit entries; signature slopes, intervals, p values, and ns as specified in Chapter 4, including TOST p 0.329 not established and permutation p 0.818; coverage 1.00, 0.73, 0.31 attached to the right species in G9; Power Map anchors 1.474 and 0.59 on the correct tiles; Savina 29 at the stated rule, 9 and 6 at 1.45, 560 at the loose rule; lifespans 4, 31, 122.5 years.

**Visual fidelity (automated where possible, else checklist):** point counts and axis ranges in G9 and G12 match their source tables; G0 renders exactly 163 tiles and 17 solid; every graph title is a conclusion sentence; the text budget of Section 5 is measured and met; species and direction palettes match Section 6 exactly; with JavaScript disabled, each chapter still shows its claim, its key numbers, and a static fallback image or table for its main graph.

**Mobile (device-tested):** every graph legible and tappable at 360 px with no horizontal scrolling anywhere; tap targets at least 44 px; no information reachable only by hover; G9 readable as stacked panels; the progress rail operable by thumb.

**Usability gate (human, launch-blocking):** Five testers who work in computational biology or hire for it, none of them briefed on the project. Four tasks, unaided, on their own phones:
1. Find the strongest aging signal in the atlas. Box: 15 seconds.
2. Explain why the human question cannot currently be answered. Box: 30 seconds after reaching Chapter 5, including naming the donor count problem.
3. On the Power Map, find a study design that meets the floor, and state its minimum detectable effect. Box: 30 seconds.
4. State whether naked mole-rat aging equivalence was proven. Correct answer: no. Box: 15 seconds after Chapter 4.
Pass bar: at least 4 of 5 testers pass each task. A failed task blocks launch; the responsible chapter is redesigned and that task retested with at least 3 fresh testers. Before any code is written, the five chapter graphs are paper-tested with printed mock graphs and the same four tasks; a story that fails on paper is fixed on paper.

**Gate:** Staging link plus a one-page build summary (Step 0 findings, which conditional exports existed, usability results, any spec deviations) goes to Reza. Public launch only on his explicit go-ahead.

## 11. Technical specification

- **Stack:** Static site, Astro preferred or an equivalent static generator settled at Step 0. One chart library. No backend, no database, no accounts. Deployable as static files under a configurable base path, working assumption `/agescape/` on the personal site, HTTPS.
- **Performance:** initial load under 1 MB excluding fallback figure images; first graph (G0) visible and legible within 2 seconds on a mid-range phone on 4G; charts render from static JSON shipped with the page; total page weight including all chapter JSON under 6 MB.
- **Accessibility:** every graph has a text twin (its expander statistics plus a data table in the provenance drawer) so no information exists only as pixels; colour is never the only encoder (fill vs hollow, solid vs dashed, and direct labels back every colour distinction); AA contrast; reduced-motion respected; full keyboard operation of rail, panel, expanders, and Power Map (tiles reachable and selectable by keyboard, readout updates identically).
- **Analytics and privacy:** optional cookieless aggregate analytics only, page views and chapter reach and Power Map tile taps in aggregate; if it cannot be cookieless, ship none. No personal data is collected. No forms in v1.

## 12. Build order and effort estimate

Estimates are planning figures for a coding agent under Reza's review, not commitments. The graph-first design removes v1's page sprawl and adds export work; net effort is comparable to v1, spent in better places.

| Step | Work | Estimate |
|---|---|---|
| 0 | Plot data audit, DATA_CONTRACT, conditional export verdicts, stack and licence decisions | 3 to 5 hours |
| 0.5 | Paper prototype of the five chapter graphs, four-task paper test, story fixes | 3 to 5 hours |
| 1 | Export script plus fidelity check, passing on real data | 4 to 6 hours |
| 2 | Page shell, visual system, Hero and Chapters 1 and 2 | 6 to 8 hours |
| 3 | Chapters 3 and 4, unit panel | 5 to 7 hours |
| 4 | Chapter 5 Power Map, Methods overlay, footer | 4 to 6 hours |
| 5 | Accessibility, mobile, and performance passes; staging deploy | 3 to 5 hours |
| 6 | Usability gate, fixes, retest, launch summary to Reza | 3 to 6 hours |

Total estimate: roughly 31 to 48 agent hours plus Reza's review at Steps 0, 0.5, and 6. Step 0.5 is not optional padding; it is the cheapest point in the project to discover the story does not land.

## 13. Decisions only Reza can make

1. Whether conditional exports C1 to C3 exist or may be produced, and, where producing one needs an analysis re-run rather than a transcription, whether that re-run is authorised at all. Default if unanswered: ship without them.
2. Canonical URL and host (personal site subpath recommended), and the final public title.
3. Code and content licences (working assumptions: MIT for code, CC BY 4.0 for report and figures).
4. Whether the Gladyshev and AP-1 validation lands before launch (adding one supporting graph to Chapter 3) or as v1.1. The chapter structure absorbs either without redesign.
5. Cookieless analytics: yes or no.
6. Public launch go or no-go, after the usability gate and his click-through.

## 14. Deferred (v2 and beyond)

Labelled interpolation on the Power Map with a methods note; UMAP explorer if embeddings are ever exported as app data; AnAge species browser with honest availability flags; Tabula Sapiens v2 or SenNet restricted-path updates, only on new committed results; an archived citable release with a DOI. None of these may be promised anywhere in the v1 app.
