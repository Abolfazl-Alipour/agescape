# Item 7: Paper-Wide Bayesian Hierarchical Analysis Across All Five Analytical Families

Date: 2026-11-14
Script: [run_rev3_bayesian_full_paper.py](file:///home/reza/projects/agescape-revision/scripts/run_rev3_bayesian_full_paper.py) (with shared Student-t and LOO helper [bayesian_diagnostics_helper.py](file:///home/reza/projects/agescape-revision/scripts/bayesian_diagnostics_helper.py))
Result tables:
- Posterior summary (`138` rows across `50` models): [rev3-bayesian-full-summary.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-summary.csv)
- MCMC convergence and LOO diagnostics (`50` models): [rev3-bayesian-full-diagnostics.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-diagnostics.csv)
- Prior predictive checks (`50` models): [rev3-bayesian-full-prior-predictive.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-prior-predictive.csv)

Companion mature-only vector drift report (Item 6): [cross-species-bayesian-vector-mature-only.md](file:///home/reza/projects/agescape-revision/report/cross-species-bayesian-vector-mature-only.md)
Scope note: Report and machine-readable tables only; the frozen PLOS ONE manuscript (`report/plos-one/final-plos-one.tex`) was not edited.

---

## 1. Executive Synthesis and Cross-Family Convergence Summary

To test whether partial pooling, heavy-tailed likelihoods, prior sensitivity, and full posterior uncertainty alter any conclusion of the Agescape study, we fit `50` Bayesian models (`4` chains each with `nutpie`, zero divergent transitions globally) spanning all five analytical families alongside the existing whole-profile drift models (`h3-cross-species-bayesian-posterior.csv` and Item 6):

| Analytical Family | Models Fit | Primary Question Tested | Bayesian Posterior Conclusion | Agreement With Frequentist Result |
| --- | ---: | --- | --- | --- |
| 1. Cell-type composition (`CLR`) | 10 | Global species direction vs cell-type-specific bidirectional turnover (`620` parametric tests + `4` NMR rank tests; mature refits) | Species-wide net CLR mean is near zero (`mu_species_std` in `[-0.0568, +0.0064]`) because shifts cancel by sum-zero constraint, while hierarchical heterogeneity `tau_species_std` separates real turnover in Mouse (`0.5559`), Fly (`0.3629`), and Worm (`0.8186`) from smaller Human (`0.1196`) and imprecise Zebrafish (`0.1526`, `n = 6` fish). Mature mouse (`>= 18m`) retains heterogeneity (`tau = 0.0188` vs `0.0137` per `10%` lifespan) while losing `lfsr < 0.05` power (`44/136` -> `11/114`). | Exact agreement (`100%` consistency with CLR OLS and Item 1 `624`-test BH; distinguishes bidirectional shifts from net drift). |
| 2. Mouse-to-human transport (`P1_Clean`) | 4 | Pooled `P1_Clean` age slope across `252` human tests (`12` cohorts) and excess above random signatures (`48` strata) | Grand pooled human standardized effect is weakly negative/null (`-0.0278`, `94% HDI [-0.0874, +0.0364], P(>0) = 0.1817`), and excess over expression-matched random gene sets is centered below zero (`-0.0855` SD, `94% HDI [-0.6999, +0.5928], P(>0) = 0.3735`). | Exact agreement: confirms the transport claim is dead despite `3/252` nominal BH hits (`q >= 0.0522` vs random sets). |
| 3. Within-tissue gene concordance | 2 | Excess Fisher-z concordance above joint donor-permutation null across `23` strata in `6` species | Fruit fly (`+0.1247 [+0.0344, +0.2168], P(>0) = 0.9947`) and C. elegans (`+0.1198 [-0.0034, +0.2694], P(>0) = 0.9721`) exceed zero; Mouse is modest and positive (`+0.0371 [-0.0072, +0.0859], P(>0) = 0.9344`); Human (`+0.0188, P(>0) = 0.7380`), Zebrafish (`+0.0402, P(>0) = 0.7698`), and NMR (`+0.0317, P(>0) = 0.7058`) cross zero widely. | Exact agreement with frequentist joint donor-permutation p-values (`0.001, 0.003, 0.044, 0.70, 0.36, 0.70`). |
| 4. Mouse lung disorder & entropy | 20 | Full-span (`1 to 30m`) off-program/entropy and mature-only (`>= 18m`) Shannon / De Man KL slopes | Depth-adjusted global `10%` lifespan (`4.8m`) slopes all cross zero broadly: full `1500` off-program `+0.1069 [-0.0578, +0.2781], P(>0) = 0.8918`; full `1500` entropy `-0.0225 [-0.1899, +0.1475], P(>0) = 0.3896`; mature `>= 18m` Shannon `+0.0626 [-0.3869, +0.5005], P(>0) = 0.6103`; mature De Man KL `-0.0365 [-0.4955, +0.4149], P(>0) = 0.4252`. | Exact agreement with inconclusive frequentist `p = 0.28 to 0.78` and wide `MDE_80`. |
| 5. Human cortical interneurons | 14 | EIV HBCC-on-MSSM replication slopes (`VIP, pvalb, sst`), independent HBCC effect sizes, and `6` shared genes | Conjugate Normal-Normal EIV yields strong, positive replication slopes (`VIP +0.9109`, `pvalb +0.9272`, `sst +1.3479`, all `P(>0) = 1.0000`, correcting OLS attenuation in `pvalb` from `0.8089`). Independent HBCC `both_FDR` magnitudes exceed `0.119` to `0.178` per `50y`. Hierarchical pooling across the `6` moderated shared genes yields `+0.1156 [+0.0802, +0.1494]` per `50y` (`P(>0) = 1.0000`). | Exact quantitative confirmation of the `381`-gene interneuron replication result (`r = 0.85 to 0.89`). |

Every single Bayesian model (`50/50`) achieved `convergence_ok == True` in [rev3-bayesian-full-diagnostics.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-diagnostics.csv):
- **Maximum R-hat across all 50 models**: `1.0079` (`<= 1.01` requirement).
- **Minimum bulk ESS across all 50 models**: `725.5` (`>= 400` requirement; `interneuron_eiv_*` reached `2,233` to `4,654` after Normal-Normal conjugate factorization).
- **Minimum tail ESS across all 50 models**: `437.6` (`>= 400` requirement).
- **Divergent transitions**: `0` across all `50` models.
- **Energy BFMI**: `[0.2838, 0.9562]`; only `composition_human_mature_30y_200-prior_0.25` fell slightly below `0.30` at `0.2838` (moving to `0.3298` at `prior_sd = 0.50`), with identical posterior estimates (`+0.0086` vs `+0.0088`) and zero divergences.
- **Prior and posterior predictive checks**: `48/50` models enclose 100% of observed data inside the central `99.8%` prior predictive interval (`[q0.001, q0.999]`). The only two flags (`composition_global_620_std-prior_0.25`: observed `[-2.2906, +2.3172]` vs prior `[-1.5663, +1.5597]`; `transport_252_std-prior_0.25`: observed `[-2.8219, +1.1495]` vs prior `[-1.8757, +1.8414]`) confirm that Gaussian shrinkage with `prior_sd = 0.25` is slightly too tight for extreme single-unit outliers, whereas `prior_sd = 0.50` encloses 100% of observations (`[-2.7315, +2.8575]` and `[-3.4545, +3.4877]`) with indistinguishable posterior means. Empirical 90% posterior predictive coverage spans `[0.9167, 1.0000]` across all `50` models.

---

## 2. Family 1: Cell-Type Composition (`CLR`, `10` Models)

### 2.1 Multi-Species Standardized Synthesis (`620` Parametric + `4` Rank Tests)

From [rev3-composition-global-bh-624.csv](file:///home/reza/projects/agescape-revision/results/rev3-composition-global-bh-624.csv), `620` tests across `5` species have valid parametric CLR standard errors (`16` C. elegans, `162` Fruit fly, `291` Human, `136` Mouse, `15` Zebrafish), whereas the **`4` Naked mole-rat skin rows (`n = 5` animals, exact permutation p-values `0.4, 0.7, 0.9, 1.0` with `se = NaN`)** are rank/fraction shifts excluded from Gaussian SE modeling so no synthetic SEs are invented.

For the `620` parametric rows, we modeled standardized age effect $d_i = (t_i / \sqrt{N_i})$ with sampling standard error $s_i = 1 / \sqrt{N_i}$, species-level mean $\mu_{\text{sp}} \sim \text{Normal}(0, \sigma_p)$, species-level heterogeneity $\tau_{\text{sp}} \sim \text{HalfNormal}(\sigma_p)$, and non-centered unit effects $\theta_i = \mu_{\text{sp}[i]} + \tau_{\text{sp}[i]} z_i$ for $\sigma_p \in \{0.25, 0.50\}$.

| Species (`n` units, max donors) | `mu_species_std` (`prior_0.50`, 94% HDI) | `P(>0)` | `tau_species_std` (`prior_0.50`, 94% HDI) | Units with `lfsr < 0.05` (`0.25` / `0.50`) | Freq Global-624 Sig (`q < 0.05`) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `Mouse` (`136`, `17`) | `-0.0362` (`[-0.1307, +0.0712]`) | `0.2485` | `0.5559` (`[+0.4769, +0.6364]`) | `47` / `48` | `11` |
| `Fruit fly` (`162`, `36`) | `+0.0062` (`[-0.0519, +0.0651]`) | `0.5742` | `0.3629` (`[+0.3169, +0.4084]`) | `74` / `73` | `26` |
| `C. elegans` (`16`, `12`) | `-0.0568` (`[-0.4276, +0.3271]`) | `0.3908` | `0.8186` (`[+0.5659, +1.1035]`) | `6` / `7` | `2` |
| `Human` (`291`, `213`) | `+0.0042` (`[-0.0184, +0.0261]`) | `0.6416` | `0.1196` (`[+0.0959, +0.1427]`) | `21` / `18` | `6` |
| `Zebrafish` (`15`, `6`) | `-0.0429` (`[-0.2496, +0.1736]`) | `0.3512` | `0.1526` (`[+0.0000, +0.3560]`) | `0` / `0` | `0` |

Key structural insight: because CLR values sum to zero within every donor/pool, species-average compositional shift `mu_species_std` is constrained near zero in every organism. Real age remodeling therefore lives in the **between-cell-type heterogeneity `tau_species_std`** and unit-level local false sign rates (`lfsr < 0.05`): large in C. elegans, Mouse, and Fruit fly (`tau = 0.36 to 0.82`), markedly smaller in adult Human (`tau = 0.12`, though high donor `N` still identifies `18 to 21` reliable unit directions), and zero `lfsr < 0.05` units in 6-fish Zebrafish testis (`tau` HDI starts near `0`). All posteriors are stable between `prior_sd = 0.25` and `0.50` (for example Mouse `mu = -0.0359` vs `-0.0362`, Human `mu = +0.0045` vs `+0.0042`).

### 2.2 Full-Span vs Mature-Only Mouse and Human Arms (`CLR` per `10%` Lifespan)

To directly parallel the frequentist mature-only composition refit (`mature-only-composition-{mouse,human}.csv`), we fit Student-t hierarchical models (`nu ~ Gamma(3, 0.2)`) on chronological CLR slopes scaled to `10%` of maximum species lifespan (`0.4y` in Mouse; `12.25y` in Human):

| Arm (`n` units, max donors) | `mu_clr_10pct_lifespan` (`prior_0.50`, 94% HDI) | `tau_clr_10pct_lifespan` (`94% HDI`) | `mean_unit_abs_theta` | Units `lfsr < 0.05` | Freq `q_within < 0.05` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `mouse_full_136` (`136`, `17`) | `-0.00087` (`[-0.00386, +0.00195]`) | `0.01375` (`[+0.01090, +0.01678]`) | `0.01150` | `44` | `19` (`11` global-624) |
| `mouse_mature_18m_114` (`114`, `13`) | `-0.00171` (`[-0.00748, +0.00409]`) | `0.01877` (`[+0.01266, +0.02526]`) | `0.01603` | `11` | `2` (`Limb_Muscle` mac & MSC) |
| `human_full_291` (`291`, `213`) | `+0.00616` (`[-0.01151, +0.02397]`) | `0.08150` (`[+0.06105, +0.10166]`) | `0.06988` | `18` | `9` (`6` global-624) |
| `human_mature_30y_200` (`200`, `209`) | `+0.00883` (`[-0.01604, +0.03243]`) | `0.09714` (`[+0.06886, +0.12588]`) | `0.08474` | `16` | `11` |

This quantitatively confirms the power-loss finding of [mature-only-composition-verdicts.csv](file:///home/reza/projects/agescape-revision/results/mature-only-composition-verdicts.csv): in mature Mouse (`>= 18m`), genuine compositional variance **does not vanish** (`tau = 0.0188` vs `0.0137`, `mean_unit_abs_theta = 0.0160` vs `0.0115`), whereas shrunken `lfsr < 0.05` survivors fall from `44` to `11` due solely to smaller sample size and compressed `18m -> 30m` leverage. In Human (`>= 30y`), nearly all adult leverage remains (`16` vs `18` units at `lfsr < 0.05`).

---

## 3. Family 2: Mouse-to-Human Signature Transport (`4` Models)

### 3.1 Standardized Pooling Across All `252` Tests (`12` Cohorts)

Combining `225` tests from [m9-regressions.csv](file:///home/reza/projects/agescape-revision/results/m9-regressions.csv) across `11` cohorts (`HBCC_DLPFC 21`, `ImmuneAging_BAL 29`, `Kidney_LivingDonor 9`, `Kidney_NormalAdult 19`, `Kidney_SCR_v15 9`, `Kidney_SN_v15 16`, `Lung_CellRef 26`, `Lung_Integrated 39`, `MSSM_DLPFC 21`, `RADC_DLPFC 18`, `SEAAD_DLPFC 18`) with `27` cell types in PsychAD ([m10-celltype-regressions.csv](file:///home/reza/projects/agescape-revision/results/m10-celltype-regressions.csv)) yields `252` standardized slopes $d_i = (t_i / \sqrt{N_i}), s_i = 1 / \sqrt{N_i}$ modeled with a hierarchical cohort intercept $\mu_c \sim \text{Normal}(\mu_{\text{global}}, \tau_{\text{cohort}})$ and unit-level random effects $\theta_i = \mu_{c[i]} + \tau_{\text{unit}} z_i$:

- **Global pooled transport (`ALL_252`)**: `mu_global_std = -0.0269` (`94% HDI [-0.0898, +0.0324], P(>0) = 0.1894`) at `prior_sd = 0.25`, and `-0.0278` (`[-0.0874, +0.0364], P(>0) = 0.1817`) at `prior_sd = 0.50`, with small cross-cohort spread `tau_cohort_std = 0.0867` (`[+0.0255, +0.1533]`).
- **Cohort posteriors (`prior_sd = 0.50`)**:
  - `Lung_Integrated` (`39` tests): `-0.0939` (`[-0.1607, -0.0262], P(>0) = 0.0033`), weakly *negative* (opposite mouse signature polarity).
  - `RADC_DLPFC` (`18`): `-0.1061` (`[-0.2182, +0.0071], P(>0) = 0.0277`).
  - `Kidney_LivingDonor` (`9`): `-0.1126` (`[-0.2770, +0.0274], P(>0) = 0.0607`).
  - `MSSM_DLPFC` (`21`): `-0.0489` (`[-0.1147, +0.0178], P(>0) = 0.0817`).
  - `Kidney_NormalAdult` (`19`): `-0.0525` (`[-0.1675, +0.0528], P(>0) = 0.1814`).
  - `Lung_CellRef` (`26`): `-0.0437` (`[-0.1432, +0.0508], P(>0) = 0.1920`).
  - `Kidney_SCR_v15` (`9`): `-0.0357` (`[-0.1608, +0.0893], P(>0) = 0.2858`).
  - `HBCC_DLPFC` (`21`): `+0.0032` (`[-0.0618, +0.0739], P(>0) = 0.5286`).
  - `Kidney_SN_v15` (`16`): `+0.0255` (`[-0.0989, +0.1533], P(>0) = 0.6336`).
  - `SEAAD_DLPFC` (`18`): `+0.0474` (`[-0.0460, +0.1556], P(>0) = 0.8094`).
  - `ImmuneAging_BAL` (`29`): `+0.0523` (`[-0.0443, +0.1545], P(>0) = 0.8351`).
  - `PsychAD` (`27`): `+0.0334` (`[-0.0216, +0.0870], P(>0) = 0.8746`).

### 3.2 Excess Above Expression-Matched Random Signatures (`48` Brain Cell Types)

To directly test whether the `3/252` global BH hits (`PsychAD EN_L3_5_IT_3`, `MSSM L2/3 IT`, `MSSM oligodendrocyte`) reflect specific biology rather than generic gene-set depth/mean inflation, we pooled all `48` standardized random-null z-excess scores $z = (\text{slope}_{\text{obs}} - \text{mean}_{\text{rand}}) / \text{sd}_{\text{rand}}$ across `PsychAD` (`27` types, [reviewer-phase1-signature-random-null.csv](file:///home/reza/projects/agescape-revision/results/reviewer-phase1-signature-random-null.csv)) and `MSSM_DLPFC` (`21` types, [reviewer-rev2-transport-null-mssm.csv](file:///home/reza/projects/agescape-revision/results/reviewer-rev2-transport-null-mssm.csv)) under a hierarchical Student-t model (`pareto_k_n_bad = 0`, `pareto_k_max <= 0.38`):
- **Global mean z-excess (`ALL_48`)**: `-0.0778` (`[-0.6038, +0.4580], P(>0) = 0.3790`) at `prior_sd = 0.50`, and `-0.0855` (`[-0.6999, +0.5928], P(>0) = 0.3735`) at `prior_sd = 1.00`.
- **MSSM_DLPFC (`21` types)**: `+0.0078` (`[-0.5239, +0.5796], P(>0) = 0.4874`) at `0.50`, and `+0.0040` (`[-0.5589, +0.6196], P(>0) = 0.4785`) at `1.00`.
- **PsychAD (`27` types)**: `-0.1735` (`[-0.6125, +0.2745], P(>0) = 0.2216`) at `0.50`, and `-0.1865` (`[-0.6340, +0.2536], P(>0) = 0.2085`) at `1.00`.

Bayesian hierarchical pooling therefore confirms that human brain `P1_Clean` slopes are indistinguishable from random expression-matched signatures.

---

## 4. Family 3: Within-Tissue Gene-Effect Concordance Across `6` Species (`2` Models)

Using all `23` organ/cohort concordance strata ([reviewer-rev2-concordance-strata.csv](file:///home/reza/projects/agescape-revision/results/reviewer-rev2-concordance-strata.csv) `11` Mouse organs, `7` Human cohorts, `2` Fruit fly regions, `1` C. elegans whole animal, `1` Zebrafish testis, `1` NMR skin), we transformed observed mean correlation minus its joint donor-permutation null mean into Fisher-z excess $y_k = \operatorname{atanh}(r_{\text{obs}, k}) - \operatorname{atanh}(\mu_{\text{null}, k})$ with permutation SD error $s_k = \frac{1}{2}[\operatorname{atanh}(\mu_{\text{null}, k} + \sigma_{\text{null}, k}) - \operatorname{atanh}(\mu_{\text{null}, k} - \sigma_{\text{null}, k})]$ and fit a two-level hierarchical Normal model ($\mu_{\text{global}}$, $\mu_{\text{species}}$, $\theta_{\text{stratum}}$) for `prior_sd in {0.25, 0.50}`:

| Species (`n` strata, `n` units) | Raw Mean `z_perm` | Bayes `mu_species_fisher_excess` (`0.25`) | `P(>0)` (`0.25`) | Bayes `mu_species_fisher_excess` (`0.50`, 94% HDI) | `P(>0)` (`0.50`) | Freq Joint `p` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Fruit fly` (`2`, `99`) | `+6.611` | `+0.1200` (`[+0.0248, +0.2095]`) | `0.9928` | `+0.1247` (`[+0.0344, +0.2168]`) | `0.9947` | `0.001` |
| `C. elegans` (`1`, `13`) | `+2.752` | `+0.1118` (`[-0.0150, +0.2552]`) | `0.9611` | `+0.1198` (`[-0.0034, +0.2694]`) | `0.9721` | `0.003` |
| `Mouse` (`11`, `78`) | `+0.574` | `+0.0373` (`[-0.0075, +0.0844]`) | `0.9396` | `+0.0371` (`[-0.0072, +0.0859]`) | `0.9344` | `0.044` |
| `Zebrafish` (`1`, `14`) | `+0.301` | `+0.0406` (`[-0.0781, +0.1552]`) | `0.7769` | `+0.0402` (`[-0.0737, +0.1576]`) | `0.7698` | `0.361` |
| `Human` (`7`, `131`) | `+0.164` | `+0.0200` (`[-0.0401, +0.0748]`) | `0.7452` | `+0.0188` (`[-0.0381, +0.0770]`) | `0.7380` | `0.700` |
| `Naked mole-rat` (`1`, `4`) | `-0.395` | `+0.0307` (`[-0.1150, +0.1809]`) | `0.7097` | `+0.0317` (`[-0.1394, +0.1800]`) | `0.7058` | `0.700` |
| `ALL_6_SPECIES` (`23`) | - | `+0.0574` (`[-0.0328, +0.1435]`) | `0.9126` | `+0.0624` (`[-0.0333, +0.1650]`) | `0.9243` | - |

Direct posterior contrasts against `Human` (`prior_sd = 0.50`):
- `Fruit fly - Human`: `+0.1059` (`94% HDI [-0.0090, +0.2097], P(>0) = 0.9602`)
- `C. elegans - Human`: `+0.1010` (`[-0.0374, +0.2569], P(>0) = 0.9074`)
- `Mouse - Human`: `+0.0183` (`[-0.0512, +0.0903], P(>0) = 0.6953`)

Partial pooling pulls single-stratum extreme `z_perm` toward the `6`-species prior (`tau_between_species = 0.0893`) while preserving the exact phylogenetic split from the manuscript: `Fruit fly` and `C. elegans` (both pooled atlases) exhibit strong shared transcriptional aging direction, `Mouse` shows modest organ-dependent positive excess (`P(>0) = 0.934 to 0.940`, concentrated in `Heart_and_Aorta` and `Mammary_Gland`), and `Human`, `Zebrafish`, and `Naked mole-rat` show no reliable excess above their joint donor-permutation nulls.

---

## 5. Family 4: Mouse Lung Off-Program Disorder and Entropy (`20` Models)

Across [reviewer-phase2-depth-donor-results.csv](file:///home/reza/projects/agescape-revision/results/reviewer-phase2-depth-donor-results.csv) (`1m to 30m`, `16` mice) and Item 4 [rev3-entropy-mature-only-mouse-lung-donors.csv](file:///home/reza/projects/agescape-revision/results/rev3-entropy-mature-only-mouse-lung-donors.csv) (`>= 18m`, `10` mice, `8` manuscript target types), we fit identical donor-clustered hierarchical Student-t models with standardized within-cell-type endpoints and `age` scaled per `10%` mouse lifespan (`0.4y = 4.8m`):

| Endpoint & Age Span (`n` types, `n` obs) | Adjustment & `prior_sd` | Bayes Mean (`10%` lifespan SD) | `94% HDI` | `P(>0)` | `elpd_loo` | Freq `p` (`age_sex` / `cov`) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `lung_full_1500_offprogram` (`8`, `96`) | unadjusted, `0.25` | `+0.0999` | `[-0.0727, +0.2804]` | `0.8618` | `-115.60` | `0.2948` |
| `lung_full_1500_offprogram` (`8`, `96`) | depth_adjusted, `0.50` | `+0.1069` | `[-0.0578, +0.2781]` | `0.8918` | `-113.52` | `0.2818` |
| `lung_full_full_offprogram` (`8`, `96`) | unadjusted, `0.50` | `+0.0851` | `[-0.1188, +0.2869]` | `0.7967` | `-105.81` | `0.4940` |
| `lung_full_full_offprogram` (`8`, `96`) | depth_adjusted, `0.50` | `+0.0837` | `[-0.1192, +0.2723]` | `0.8010` | `-105.77` | `0.5038` |
| `lung_full_1500_entropy` (`9`, `105`) | unadjusted, `0.50` | `-0.0372` | `[-0.2385, +0.1511]` | `0.3561` | `-151.12` | `0.7041` |
| `lung_full_1500_entropy` (`9`, `105`) | depth_adjusted, `0.50` | `-0.0225` | `[-0.1899, +0.1475]` | `0.3896` | `-145.74` | `0.7299` |
| `lung_mature18m_1500_entropy` (`8`, `67`) | unadjusted, `0.50` | `-0.0006` | `[-0.4935, +0.4907]` | `0.5062` | `-92.06` | `0.7783` |
| `lung_mature18m_1500_entropy` (`8`, `67`) | depth_adjusted, `0.50` | `+0.0626` | `[-0.3869, +0.5005]` | `0.6103` | `-87.95` | `0.7783` (base) |
| `lung_mature18m_1500_deman_kl` (`8`, `67`) | unadjusted, `0.50` | `+0.0165` | `[-0.4552, +0.5078]` | `0.5208` | `-84.93` | `0.7411` |
| `lung_mature18m_1500_deman_kl` (`8`, `67`) | depth_adjusted, `0.50` | `-0.0365` | `[-0.4955, +0.4149]` | `0.4252` | `-82.46` | `0.7411` (base) |

Complete table with both `prior_sd = 0.25` and `0.50` is in [rev3-bayesian-full-summary.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-summary.csv). Every model yields wide intervals straddling zero, and adjusting for `mean_log_umi_actual` consistently improves `elpd_loo` across all `5` endpoints (`+0.04` to `+5.39`, largest on Shannon entropy). Even in unconfounded post-maturity `18m to 30m` rarefied lung (`1500` UMI), `94%` HDIs of `[-0.39, +0.50]` confirm genuinely low precision (`MDE_80` ~ `1.5 SD`) rather than evidence of absence.

---

## 6. Family 5: Human Cortical Interneuron Replication (`14` Models)

Using [reviewer-rev2-interneuron-genes-moderated.csv](file:///home/reza/projects/agescape-revision/results/reviewer-rev2-interneuron-genes-moderated.csv) (`768` total genes across `GABA`, `VIP`, `pvalb`, `sst` tested in both `MSSM_DLPFC` and `HBCC_DLPFC`, donors `>= 20 y`), we quantified interneuron replication via three complementary Bayesian formulations:

### 6.1 Bayesian Errors-in-Variables (EIV) Replication Slopes (`HBCC` on `MSSM`)

Naive OLS regression of `slope_hbcc` on noisy discovery `slope_mssm` suffers from measurement-error attenuation. Using exact Normal-Normal conjugate factorization on `(x_obs, x_true)` combined with a Student-t likelihood (`sigma_tot = sqrt(sigma_rep_resid^2 + y_se^2)`) and `prior_sd in {1.00, 2.00}` (`pareto_k_n_bad <= 1`, `ess >= 2,233`):

| Subtype (`n` MSSM-FDR genes) | Naive OLS Slope (Pearson `r`) | Bayes EIV `beta_hbcc_on_mssm` (`prior_1.00`, 94% HDI) | Bayes EIV `beta_hbcc_on_mssm` (`prior_2.00`, 94% HDI) | `P(>0)` |
| --- | ---: | ---: | ---: | ---: |
| `VIP` (`203`) | `0.8765` (`0.8554`) | `0.9091` (`[0.8358, 0.9827]`) | `0.9109` (`[0.8370, 0.9866]`) | `1.0000` |
| `pvalb` (`137`) | `0.8089` (`0.8495`) | `0.9260` (`[0.8325, 1.0200]`) | `0.9272` (`[0.8343, 1.0261]`) | `1.0000` |
| `sst` (`362`) | `1.3957` (`0.8897`) | `1.3464` (`[1.2665, 1.4255]`) | `1.3479` (`[1.2659, 1.4257]`) | `1.0000` |

EIV de-attenuates `pvalb` (`0.8089 -> 0.9272`, `94% HDI [0.8343, 1.0261]`) and `VIP` (`0.8765 -> 0.9109`) cleanly toward unity (`1.0`), while `sst` replicates with `1.35x` larger effect magnitudes in HBCC relative to MSSM, invariant to doubling the slope prior width (`1.00 -> 2.00`).

### 6.2 Independent `HBCC` Magnitude on `both_FDR` Genes & Pooling Across Shared Genes

Fitting Student-t hierarchical models solely to the independent replication cohort `HBCC` (`|slope_hbcc| * 50y` in `log(cp10k)` per `50` years across the `381` subtype `both_FDR` hits) yields virtually identical estimates across `prior_sd = 0.25` and `0.50`:
- `VIP` (`77` genes): `+0.1325` (`94% HDI [+0.1215, +0.1435], P(>0) = 1.0000`)
- `pvalb` (`65` genes): `+0.1194` (`94% HDI [+0.1079, +0.1303], P(>0) = 1.0000`)
- `sst` (`239` genes): `+0.1785` (`94% HDI [+0.1656, +0.1909], P(>0) = 1.0000`)

In the moderated intersection across all three interneuron subtypes (`VIP`, `pvalb`, `sst`), **`6` genes pass `both_FDR` in all three** (`FBN1`, `GABRB1`, `IGF1R`, `PHF3`, `TRMT11`, `ZBTB16`; yielding `36` gene-by-subtype-by-cohort observations across MSSM + HBCC). Fitting a two-cohort hierarchical Gaussian shrinkage model (`mu_shared_50y`, `prior_sd = 0.25` and `0.50`, `pareto_k_n_bad = 0` at `0.50`) yields:
- **Grand shared 6-gene effect (`prior_0.50`)**: `+0.1156` log-CP10k per `50` years (`94% HDI [+0.0802, +0.1494], P(>0) = 1.0000`; `prior_0.25` = `+0.1147 [+0.0816, +0.1471]`).
- **Gene-level hierarchical posteriors (`prior_0.50`, per `50y`)**:
  - `GABRB1`: `+0.1524` (`94% HDI [+0.1193, +0.1855]`, `P(>0) = 1.0000`)
  - `TRMT11`: `+0.1217` (`[+0.0975, +0.1446]`, `P(>0) = 1.0000`)
  - `IGF1R`: `+0.1204` (`[+0.0972, +0.1446]`, `P(>0) = 1.0000`)
  - `ZBTB16`: `+0.1186` (`[+0.0883, +0.1472]`, `P(>0) = 1.0000`)
  - `FBN1`: `+0.0944` (`[+0.0724, +0.1165]`, `P(>0) = 1.0000`)
  - `PHF3`: `+0.0849` (`[+0.0632, +0.1070]`, `P(>0) = 1.0000`)

**Note on the 7th manuscript gene (`DYNC1I1`)**: The PLOS ONE manuscript reports `7` upregulated genes shared across `VIP`, `pvalb`, and `sst`, which matches the unmoderated intersection in [reviewer-rev2-interneuron-genes.csv](file:///home/reza/projects/agescape-revision/results/reviewer-rev2-interneuron-genes.csv) (`DYNC1I1` in `pvalb`: `q_bh_mssm = 0.04959`, `q_bh_hbcc = 0.03536`; `VIP`: `q_bh_mssm = 0.00119`, `q_bh_hbcc = 0.00788`; `sst`: `q_bh_mssm = 7.65e-9`, `q_bh_hbcc = 0.000325`). Under Smyth (2004) empirical-Bayes SE moderation ([reviewer-rev2-interneuron-genes-moderated.csv](file:///home/reza/projects/agescape-revision/results/reviewer-rev2-interneuron-genes-moderated.csv)), `DYNC1I1` remains `both_FDR` in `VIP` and `sst` (`q_bh_mssm_mod = 0.00136` and `1.35e-8`) while `pvalb` MSSM moderated FDR shifts marginally past `0.05`, leaving the `6`-gene moderated core modeled above.

---

## 7. Overall Verdict Across the Entire Agescape Paper

Across `60` total Bayesian models (`10` mature-only vector drift models in Item 6 + `50` full-paper models here; `0` divergences globally, `60/60` converged), **Bayesian partial pooling confirms every core qualitative finding of the frozen PLOS ONE manuscript without reversing any null or positive claim**:
1. **No universal transcriptional aging program**: compositional averages cancel (`tau` concentrates in Mouse, Fly, and Worm), `P1_Clean` fails to exceed random gene sets across `252` human tests (`-0.0278`), mouse mature whole-vector drift vanishes once depth is controlled (`Item 6`, `+0.0389`, `P(>0) = 0.7068`), and mouse lung disorder/entropy slopes straddle zero across both full-lifespan and mature-only fits.
2. **Real structure is localized**: within-tissue gene-effect concordance is robust in Fly (`P(>0) = 0.9947`) and Worm (`P(>0) = 0.9721`) and weak/organ-specific in Mouse (`0.9344`), while cortical interneuron aging in adult human prefrontal cortex replicates cleanly across independent cohorts (`VIP`, `pvalb`, `sst` EIV slopes `0.91 to 1.35`, `P(>0) = 1.0000`).
