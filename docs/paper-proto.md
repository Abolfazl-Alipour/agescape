# Paper Prototype Kit — Step 0.5 (for Reza to run with 5 testers)

**Purpose (spec §12):** test the five chapter graphs on paper with the four usability tasks BEFORE UI code, because a story that fails on paper fails harder in a browser. **This is a launch prerequisite, not optional padding.**

**How to run:** print this file (the mock blocks read fine in monospace at A4). Recruit 5 people who work in computational biology or hire for it; do not brief them on the project. Hand them the sheets in chapter order, read each task aloud, time silently, record on the score sheet. Do not answer questions during a task; note every question asked.

---

## MOCK G0/G1 — Chapter 1 ("Of 163 cell types tested, aging is powered in 17, all in lung and kidney")

```
Tissue            units (each ■ = one cell type; solid = powered, □ = ineligible)
Lung              ■ ■ ■ ■ ■ ■ ■ ■ ■ □ □ □ □ □ □ □ □ □ □ □ □     Lung 9 of 21
Kidney            ■ ■ ■ ■ ■ ■ ■ ■ □ □ □ □ □ □ □ □ □ □ □ □ □ □ □  Kidney 8 of 23
Marrow            □ □ □ □ □ □ □ □ □ □ □ □ □ □ □ □ □ □            Marrow 0 of 18
Spleen            □ □ □ □ □ □ □ □ □ □ □ □                        Spleen 0 of 12
Liver             □ □ □ □ □ □ □ □ □                              Liver 0 of 9
... (11 more rows, all hollow)

Ranking (powered units only)
Lung · adventitial cell          |██████████████████████████████████████████████████ 55
Kidney · prox. convoluted tubule |██████████████████████████████████████████████    53
Lung · fibroblast                |██████████████████████████████████                36
Lung · alveolar macrophage       |█████████████████████████████████                 35
  (13 more bars, 28 → 0)
```
*(Tile counts per tissue and per-row "X of Y" come from atlas.json at build time; this mock uses the real totals: 163 tiles, 17 solid, all in two rows.)*

## MOCK G2/G3/G4 — Chapter 2 ("The top signal is 12 times its own permutation null, and the pipeline agrees with the original paper")

```
G2 Calibration (number line)
0 ──────────────●────────────────────────────┬───── sex-effect floor: 0
               4.6                          55 ┘ observed: 55
               permutation null mean         (permutation, shuffled ages)

G3 Agreement (lollipop, rho 0 → 1)
13 agree  ▮▮▮▮▮▮▮▮▮▮▮▮▮▮▮
4 partial ▮▮▮▮
           ●           ●      ●  ●    ●                ●
0 |---------|---------|---------|---------|---------|---- rho
           0.2       0.4       0.6       0.8  0.866(Kidney PCT)

   [scatter inset] our log2FC vs paper age coefficient, 13,019 genes, rho 0.866

G4 Design grid (mice per age group; hatched = single-sex column)
                 1m    3m    18m   21m   24m   30m
female             0  ▒▒▒2▒▒▒ ▒▒▒2▒▒▒ ▓▓▓2▓▓▓     0     0
male             ▓▓2▓▓ ▓▓▓5▓▓▓ ▓▓2▓▓▓     0 ▓▓4▓▓▓ ▓▓4▓▓▓
                 ^male-only            ^female-only      ^male-only
                      └─ outlined: "the only fair fight" ─┘
```

## MOCK G6/G7 — Chapter 3 ("261 genes move with age: 167 up, 94 down")

```
Program per unit (stacked; ■ up=vermillion, ■ down=blue-grey)
Lung · adventitial          ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 34u 21d
Kidney · prox. conv. tubule ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    30u 23d
   ... (15 more bars; zero-DE units drawn as empty tracks labelled 0)

Gene × unit matrix (top 40 genes; ● = in program, colour = direction)
         Lung:adv  Lung:fib  Kid: PCT  Kid:DCT  ... (17 columns)
Ccl2        ●           ●         ●                   (shared)
Spp1        ●                     ●                   (shared)
Xist        ●                                         (private)
```

## MOCK G9/G10/G11/G12 — Chapter 4 ("In mouse, the program rises with age. In the naked mole-rat, over the span tested, it does not. Whether that difference is real is not yet settled.")

```
G9 Small multiples (shared y axis; line + 95% CI band; coverage meter under each)
Mouse (solid)      Human (dashed, LOW-POWERED watermark)   NMR
score │        ●      score │  ●                          score │ ●
      │      ●   ●          │    ●  ●                           ●
      │    ●       ●        │ ●      ●                         ●   ●
      │  ●           ●      │        (band spans everything)  │ band spans 0
      └────────────────        └─────────────                  └─────────────
      slope 0.235 [0.067, 0.403]  n=4, slope −0.21 [−1.37, 0.95]  slope 0.017 [−0.257, 0.290]
      p=0.008  n=23               coverage 0.73                   p=0.858  n=5  coverage 0.31
      coverage 1.00                                                 tested span 1.1–11.3 y

G10 The honest interval (forest plot; bold vertical line at 0; band = ±δ around NMR row)
                                  │
Mouse slope            ━━━━━●━━━━━┃━━━━━                          CI excludes 0
NMR slope        ━━●━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━                CI spans 0
                                  │ ▒▒▒ ±δ band (δ=0.059) ▒▒▒ → equivalence NOT established
Slope difference      ━━━━━●━━━━━━┿━━●━━━━                         CI just touches 0 (p=0.064)

G11 Lifespan context (years, linear)
Mouse  ■■░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 4 y   (tested to 2.5 y)
NMR    ■■■■■■■■■■■░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 31 y  (tested to 11.3 y)
Human  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■░░░░░░░░░░░░░░░ 122.5 y

G12 Savina curve (DE genes vs |log2FC| cutoff)
count │●
(log) │  ●
      │    ●        ● paper rule (|log2FC|≥1.0): 29 genes — paper reported 0
      │      ●      ● converted cutoff 1.45: 9 (p) / 6 (FDR)
      │        ● ●  ● loose rule (FDR≥0.25 @0.25): 560
      └────────────────────── cutoff
```

## MOCK G13 — Chapter 5 ("With 4 donors, only enormous effects could be seen. Here is the design that could see real ones.")

```
Minimum detectable 30-year log2FC (dark = only huge effects detectable · light = sensitive)
            donors →   4     6     8    10    12    16    20    30    50
age SD 20  │  ░░░ ░░░ ░░░ ░░░ ░░░ ░░░ ░░░ ░░░ ░░░      (row SD=20: sensitive)
age SD 16  │  ░░░ ░░░ ░▒▒ ░▒▒ ░▒▒ ░▒▒ ░▒▒ ░▒▒ ░▒▒
age SD 12  │  ▒▒▒ ▒▒▒ ▒▒▓ ▒▒▓ ░▒▒ ░▒▒ ░▒▒ ░▒▒ ░▒▒      ← eligibility floor (≥6 donors, ≥12 SD)
age SD 8   │  ▓▓▓ ▓▓▓ ▓▓▓ ▓▓▓ ▓▓▓ ▒▓▓ ▒▓▓ ▒▓▓ ▒▓▓      outlined tiles = meet floor
                ✕ Tabula Sapiens skin: 4 donors, SD 12 → 1.474 (dispersion 0.1), LOW-POWERED
Readout: [donors 6] [SD 12] [dispersion 0.1] → min detectable 30y log2FC = 0.838 → MEETS FLOOR
(3 panels: dispersion 0.1 / 0.2 / 0.4 — values from power-grid.json only)
```

---

## The four tasks (read aloud, unaided; time silently)

1. "Find the strongest aging signal in the atlas." — **Box 15 s.** Pass: names Lung adventitial cell (or points at the longest bar / a solid tile in the Lung row).
2. "Explain why the human question can't currently be answered." (after Chapter 5) — **Box 30 s.** Pass: mentions donor count (4) or points at the ✕ tile / LOW-POWERED chip.
3. "On this Power Map, find a study design that meets the floor and say what effect size it could detect." — **Box 30 s.** Pass: selects an outlined tile (≥6 donors, ≥12 SD) and reads a min-detectable value off the readout.
4. "Was naked-mole-rat aging equivalence proven?" (after Chapter 4) — **Box 15 s.** Pass: "No" (equivalence not established / CI extends beyond the band).

**Pass bar (spec §10):** ≥4 of 5 testers pass EACH task. Any failed task blocks launch; redesign the responsible chapter and retest that task with ≥3 fresh testers.

## Score sheet (copy per tester)

| Tester # | Task 1 (15s) time/pass | Task 2 (30s) time/pass | Task 3 (30s) time/pass | Task 4 (15s) time/pass | Questions asked (verbatim) | Confusion points |
|---|---|---|---|---|---|---|
| 1 | | | | | | |
| 2 | | | | | | |
| 3 | | | | | | |
| 4 | | | | | | |
| 5 | | | | | | |

**Result:** each task passed by __/5. Gate verdict: PASS / FAIL (any task <4 = FAIL → redesign chapter, retest with ≥3 fresh).

**Outcome of the run:** ______________________  **Date:** ______  **Signature:** ______
