# Web Application Technology Stack & Architecture Decision

**Milestone:** M8  
**Component:** Conserved Aging Explorer Web Application (`app/`)  
**Deployment Target:** Cloudflare Pages / Vercel Static Hosting (Static client-side architecture)

---

## 1. Stack Selection & Biological Rationale

### Choice: Modern Client-Side Vanilla Web Application with Apache Arrow / Parquet Data Layer
- **No Heavy Framework Overhead:** Vanilla ES modules with high-performance CSS and SVG charting engines provide instantaneous load times (<500ms), zero dependency build vulnerabilities, and maximal longevity without framework bitrot.
- **Data Engine:** Apache Arrow / Parquet data formats for zero-copy in-memory querying. Numeric gene expression matrices and cell coordinates are stored in compact columnar format (`float16`), well below the pre-registered budget ceilings (total assets ≤ 300 MB; per-view payload ≤ 25 MB).
- **Security & Hosting:** Deployable as an immutable static site on Cloudflare Pages / Netlify / Vercel class hosting. No backend server execution, zero server-side state, zero credentials, zero attack surface.

---

## 2. Feature Prioritization & Information Architecture

1. **View 1 (Lead View): Signature Longevity-Outlier Contrast**
   - Direct interactive comparison of conserved aging program activation across normalized lifespan for **Mouse**, **Human**, and **Naked Mole-Rat**.
   - Includes interactive TOST equivalence band toggles ($\delta = \pm 0.15$), slope standard errors, and life-stage translation slider.
2. **View 2: Gene-Level Explorer & Query**
   - Searchable gene autocomplete covering all 1:1 mammalian orthologs.
   - Cross-species trajectory charts, effect sizes ($log_2\text{FC}$), 95% bootstrap confidence intervals, and conservation classifications.
3. **View 3: Conserved vs Divergent Program Panels**
   - Interactive data table of the headline divergence artifact (`results/m5-divergence-table.csv`).
   - Deep dives into AP-1 stress activation, inflammatory SASP, ECM remodeling, and mitochondrial decline.
4. **View 4: Atlas & Tissue Explorer**
   - UMAP embeddings colored by categorical age cohorts with cell composition toggle.
   - Per-tissue demographic distributions (3m vs 18m).
5. **View 5: Caveats & Methodological Transparency**
   - Direct link from every view to the comprehensive caveats documentation: sex×age droplet imbalance, cross-sectional vs longitudinal limits, and NMR culture constraints.
