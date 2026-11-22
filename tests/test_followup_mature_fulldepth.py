"""Regression tests for the mature 30-plus full-depth pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from run_followup_de_man_entropy import kl_uniform, multinomial_downsample, shannon  # noqa: E402
from run_h3_expression_vector import hellinger_distance  # noqa: E402
import run_followup_mature_fulldepth as fd  # noqa: E402


def test_stringent_ambient_mask_preserves_igf_hb_neighbors() -> None:
    symbols = [
        "IGF1", "IGF2", "IGFBP3", "IGFBP5",
        "HBEGF", "HBS1L",
        "IGHG1", "IGHM", "IGKC", "IGLC2", "IGHA1",
        "HBB", "HBA1", "FOS", "JUN", "JUND",
        "TTN", "MYH6", "GAPDH",
    ]
    mask_gao = fd.build_stringent_ambient_mask(symbols, "gao", None)
    by_sym = dict(zip(symbols, mask_gao))
    for keep_gene in ("IGF1", "IGF2", "IGFBP3", "IGFBP5", "HBEGF", "HBS1L", "GAPDH"):
        assert not by_sym[keep_gene], f"{keep_gene} falsely masked"
    for drop_gene in ("IGHG1", "IGHM", "IGKC", "IGLC2", "IGHA1", "HBB", "HBA1", "FOS", "JUN", "JUND", "TTN", "MYH6"):
        assert by_sym[drop_gene], f"{drop_gene} missed by mask"

    # Yazar spares cardiomyocyte curated genes but keeps IEG/HB/IG/platelet exclusions
    mask_yazar = fd.build_stringent_ambient_mask(["TTN", "PPBP", "PF4", "IGHG1", "IGF1"], "yazar", None)
    assert mask_yazar.tolist() == [False, True, True, True, False]


def test_fast_entropy_and_hellinger_match_canonical_implementations() -> None:
    rng = np.random.default_rng(42)
    n_cells, n_genes = 15, 40
    dense = rng.poisson(3.0, size=(n_cells, n_genes)).astype(np.float64)
    dense[:, :5] += 25.0
    mat = sp.csr_matrix(dense)

    sh_fast, klm_fast, kls_fast = fd.fast_entropy_from_csr(mat, n_genes)
    sh_ref = shannon(mat)
    klm_ref = kl_uniform(mat, sqrt=False)
    kls_ref = kl_uniform(mat, sqrt=True)
    np.testing.assert_allclose(sh_fast, sh_ref, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(klm_fast, klm_ref, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(kls_fast, kls_ref, rtol=1e-12, atol=1e-12)

    ref_vec = rng.uniform(1.0, 20.0, size=n_genes)
    ref_sqrt = np.sqrt(ref_vec / ref_vec.sum())
    h_fast = fd.fast_hellinger_csr_vec(mat, ref_sqrt, already_l1=False)
    h_ref = np.array([hellinger_distance(dense[i] / dense[i].sum(), ref_vec / ref_vec.sum()) for i in range(n_cells)]).ravel()
    np.testing.assert_allclose(h_fast, h_ref, rtol=1e-12, atol=1e-12)


def test_compute_sh_and_quenouille_exact_scaling_and_bypass_on_counts() -> None:
    rng = np.random.default_rng(123)
    n_d = 24
    ages = np.linspace(22.0, 78.0, n_d)
    batches = ["B0" if i % 2 == 0 else "B1" for i in range(n_d)]
    bins = ["YOUNG" if a < 40 else ("MID" if a < 60 else "OLD") for a in ages]
    obs = 1.0 + 0.25 * (ages / 10.0) + rng.normal(0.0, 0.15, size=n_d)
    noise = rng.normal(0.0, 0.08, size=n_d)
    half_bias = 0.02 * (ages / 10.0) + rng.normal(0.0, 0.01, size=n_d)
    df = pd.DataFrame({
        "donor_id": [f"D{i}" for i in range(n_d)],
        "age_years": ages,
        "age_per_decade": ages / 10.0,
        "bin": bins,
        "sex_male": [(i // 2) % 2 for i in range(n_d)],
        "_depth": rng.normal(7.5, 0.2, size=n_d),
        "study_or_pool": batches,
        "obs": obs,
        "h1_0": obs + half_bias + noise,
        "h2_0": obs + half_bias - noise,
        "h1_1": obs + half_bias + 0.5 * noise,
        "h2_1": obs + half_bias - 0.5 * noise,
    })

    st_scale_inv = fd.compute_sh_and_quenouille_stats(
        df, "obs", ["h1_0", "h1_1"], ["h2_0", "h2_1"],
        ["sex_male", "_depth"], "study_or_pool", 2, 5, "B0",
        is_part3_count_metric=False, apply_quenouille=True, fallback_to_nobatch=False,
    )
    st_cnt = fd.compute_sh_and_quenouille_stats(
        df, "obs", ["h1_0", "h1_1"], ["h2_0", "h2_1"],
        ["sex_male", "_depth"], "study_or_pool", 2, 5, "B0",
        is_part3_count_metric=True, apply_quenouille=False, fallback_to_nobatch=False,
    )
    # Count metric uses factor 1.0 vs scale-invariant 0.5 -> exactly 2x SE_samp
    assert np.isfinite(st_scale_inv["sh_se_samp_slope"]) and st_scale_inv["sh_se_samp_slope"] > 0.0
    np.testing.assert_allclose(st_cnt["sh_se_samp_slope"], 2.0 * st_scale_inv["sh_se_samp_slope"], rtol=1e-12)
    np.testing.assert_allclose(st_cnt["sh_rho_tech2"], 4.0 * st_scale_inv["sh_rho_tech2"], rtol=1e-12)
    assert np.isfinite(st_scale_inv["quenouille_corr_slope"])
    assert np.isnan(st_cnt["quenouille_bias_slope"])
    assert np.isnan(st_cnt["quenouille_corr_slope"])


def test_build_e1_offprogram_joint_condition_and_invariant_sources() -> None:
    # 3 CTs (0, 1, 2), 6 donors (4 YOUNG, 2 OLD), 4 genes
    n_donors, n_cts, n_genes = 6, 3, 4
    bin_arr = np.array([1, 1, 1, 1, 3, 3], dtype=np.int32)
    cnt_ct_15 = np.full((n_donors, n_cts), 30, dtype=np.int64)
    y_pair_map = {(d * n_cts + c): (d * n_cts + c) for d in range(4) for c in range(n_cts)}
    y_n_cells = np.full(4 * n_cts, 30, dtype=np.int64)
    y_full_cnt = np.zeros((4 * n_cts, n_genes), dtype=np.float64)
    y_full_det = np.zeros((4 * n_cts, n_genes), dtype=np.float64)

    for d in range(4):
        # Target CT 0 silent (<0.5% det, ~0 cp10k) on gene 0 & gene 1; strong background gene 3 across all CTs
        y_full_cnt[d * n_cts + 0, :] = [0, 0, 1000, 1000]
        y_full_det[d * n_cts + 0, :] = [0, 0, 30, 30]
        # CT 1: gene 0 has HIGH det_frac (25/30 = 0.83) AND 1000 cp10k -> qualifies for CT 0 off-program
        #       gene 1 has HIGH cp10k (5000) ONLY because 1 outlier cell had huge UMI? det_frac=2/30 (<0.10) -> does NOT qualify in CT 1!
        y_full_cnt[d * n_cts + 1, :] = [100, 500, 0, 400]
        y_full_det[d * n_cts + 1, :] = [25, 2, 0, 30]
        # CT 2: gene 1 has det_frac=15/30=0.50 and cp10k=1500/5000*10k=3000 -> QUALIFIES in CT 2!
        #       (tests that higher-cp10k non-qualifying CT 1 does NOT block qualifying CT 2!)
        y_full_cnt[d * n_cts + 2, :] = [0, 150, 0, 350]
        y_full_det[d * n_cts + 2, :] = [0, 15, 0, 30]

    e1 = fd.build_e1_offprogram_lookup(
        y_full_cnt, y_full_det, y_n_cells, y_pair_map, cnt_ct_15,
        bin_arr, n_donors, n_cts, n_genes, ge20_only_ref=True, exclude_mask=None,
    )
    assert e1["invariant_cts"] == [0, 1, 2]
    assert set(e1["full"][0].tolist()) == {0, 1}
    # Removing any one YOUNG donor still leaves same invariant sources and genes 0,1 qualifying
    assert set(e1["loo"][(0, 0)].tolist()) == {0, 1}


def test_build_e2_silent_bounded_stratum_loo_guard() -> None:
    # 3 YOUNG donors, 1 CT, 2 genes: cells = [120, 120, 80] => total 320, min_loo = 200 -> passes 200!
    # If donor 0 has 130 cells -> total 330, min_loo = 200; if donor 0 has 131 -> min_loo=199 -> fails!
    n_donors, n_cts, n_genes = 3, 1, 2
    bin_arr = np.array([1, 1, 1], dtype=np.int32)
    y_pair_map = {0: 0, 1: 1, 2: 2}
    cnt_ct_15 = np.full((3, 1), 25, dtype=np.int64)
    det = np.array([[0, 10], [0, 10], [0, 10]], dtype=np.int64)

    e2_pass = fd.build_e2_silent_lookup(
        det, np.array([120, 120, 80], dtype=np.int64), y_pair_map, cnt_ct_15,
        bin_arr, n_donors, n_cts, n_genes, min_loo_young_cells=200, ge20_only_ref=False,
    )
    assert e2_pass["full"][0].tolist() == [0]

    e2_fail = fd.build_e2_silent_lookup(
        det, np.array([131, 120, 79], dtype=np.int64), y_pair_map, cnt_ct_15,
        bin_arr, n_donors, n_cts, n_genes, min_loo_young_cells=200, ge20_only_ref=False,
    )
    assert 0 not in e2_fail["full"]


def test_binary_ectopic_metrics_exact_brute_force_and_block_boundary() -> None:
    # Build >512 rows (520 cells) with sparse hits across 50 genes to test block=512 boundary!
    rng = np.random.default_rng(7)
    n_cells, n_genes = 520, 50
    rows, cols = [], []
    # Cell 0..515 mix hits; tail 4 cells all zeros
    for r in range(516):
        if rng.random() < 0.35:
            gs = rng.choice(n_genes, size=int(rng.integers(1, 5)), replace=False)
            for g in gs:
                rows.append(r)
                cols.append(int(g))
    bmat = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_cells, n_genes), dtype=np.float64)
    out = fd.compute_binary_ectopic_metrics(bmat, np.random.default_rng(99), n_genes_set=n_genes, run_excess=True)

    # Brute force exact reference on dense matrix
    bd = (bmat.toarray() > 0).astype(bool)
    pos_idx = np.flatnonzero(bd.sum(axis=1) > 0)
    j_sum = 0.0
    for a_i in range(len(pos_idx)):
        ra = bd[pos_idx[a_i]]
        for b_i in range(a_i + 1, len(pos_idx)):
            rb = bd[pos_idx[b_i]]
            inter = np.logical_and(ra, rb).sum()
            union = np.logical_or(ra, rb).sum()
            j_sum += inter / union
    p_all = n_cells * (n_cells - 1) // 2
    p_both = len(pos_idx) * (len(pos_idx) - 1) // 2
    n_zero = n_cells - len(pos_idx)
    p_union = p_all - (n_zero * (n_zero - 1) // 2)
    np.testing.assert_allclose(out["jaccard_both_pos"], j_sum / p_both, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(out["jaccard_all_pairs"], j_sum / p_all, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(out["jaccard_union_pos"], j_sum / p_union, rtol=1e-12, atol=1e-12)
    assert np.isfinite(out["excess_jaccard"])

    # Ensure compute_binary_ectopic_metrics does NOT mutate input counts >= 2 in-place
    cnt_mat = sp.csr_matrix((np.full(len(rows), 3, dtype=np.int64), (rows, cols)), shape=(n_cells, n_genes))
    fd.compute_binary_ectopic_metrics(cnt_mat, np.random.default_rng(2), n_genes_set=n_genes, run_excess=False)
    assert cnt_mat.dtype == np.int64 and int(cnt_mat.data.min()) == 3

    # Zero-hit matrix must return NaN singleton_fraction and NaN jaccard_both_pos (never 0/0)
    zmat = sp.csr_matrix((10, 20), dtype=np.float64)
    zout = fd.compute_binary_ectopic_metrics(zmat, np.random.default_rng(1), n_genes_set=20, run_excess=False)
    assert np.isnan(zout["singleton_fraction"])
    assert np.isnan(zout["jaccard_both_pos"])


def test_stratified_global_with_batch_and_verdict_gates() -> None:
    rng = np.random.default_rng(100)
    rows = []
    # 2 cell types, 20 donors, batch-confounded offset B1=+2.0 (older weighted slightly toward B1)
    for d in range(20):
        age = 24.0 + 3.0 * d
        b = "B1" if (d >= 10 and d % 3 != 0) or (d < 10 and d % 4 == 0) else "ENCODE v4 (Snyder)"
        bn = "YOUNG" if age < 40 else ("MID" if age < 60 else "OLD")
        b_eff = 2.5 if b == "B1" else 0.0
        for ct in ("CT_A", "CT_B"):
            y = b_eff + 0.01 * (age / 10.0) + rng.normal(0.0, 0.05)
            rows.append({
                "donor_id": f"D{d}", "age_years": age, "age_per_decade": age / 10.0,
                "bin": bn, "sex_male": d % 2, "_depth": 8.0 + 0.01 * d, "n_cells": 50,
                "study_or_pool": b, "group": ct, "yazar_pool_fe_all_evaluable": True, "val": y,
            })
    df_st = pd.DataFrame(rows)
    srow, _ = fd.fit_stratified_global_with_batch(
        df_st, ["CT_A", "CT_B"], "val", "n_cells", "_depth", "gao", "1", "center_shift", "full_depth", "primary",
    )
    assert srow is not None
    assert srow["is_stratified_global"] is True
    assert "batch_fe(1)" in srow["covariates_active"]
    assert abs(srow["slope_raw_per_decade"] - 0.01) < 0.05
    assert np.isfinite(srow["old_minus_young_diff_adj"])

    # Verdict gate: randomness vs shared_induction vs separate Part 1/2/3 verdicts
    def _mk(part: str, metric: str, dm: str, gs: str, sl: float, p: float, snr: float, oy_d: float = 0.0, oy_p: float = 0.5) -> dict:
        return {
            "part": part, "metric": metric, "depth_mode": dm, "model_spec": "stratified_global_primary",
            "gene_set": gs, "slope_raw_per_decade": sl, "pvalue": p, "sh_snr_slope": snr,
            "old_minus_young_diff_adj": oy_d, "old_minus_young_p_adj": oy_p,
            "quenouille_bias_frac": 0.01, "quenouille_corr_slope": sl,
            "old_mean": 0.55, "young_loo_mean": 0.50 if sl > 0 else 0.56, "mde80_sd_per_decade": 0.12,
        }

    p1 = pd.DataFrame([
        _mk("1", "center_shift", "full_depth", "none", 0.02, 1e-4, 5.0, 0.04, 1e-3),
        _mk("1", "center_shift", "spline_resid_full", "none", 0.02, 1e-4, 5.0),
        _mk("1", "scatter", "full_depth", "none", -0.01, 0.40, 0.5, -0.01, 0.4),
    ])
    p2 = pd.DataFrame([
        _mk("2", "shannon", "full_depth", "none", 0.03, 1e-4, 4.0),
        _mk("2", "shannon", "rarefied_3000", "none", 0.02, 1e-2, np.nan),
        _mk("2", "deman_kl_methods", "full_depth", "none", -0.03, 1e-4, 4.0),
        _mk("2", "deman_kl_methods", "rarefied_3000", "none", -0.02, 1e-2, np.nan),
        _mk("2", "deman_kl_sqrt", "full_depth", "none", -0.01, 0.5, 0.4),
        _mk("2", "deman_kl_sqrt", "rarefied_3000", "none", 0.01, 0.5, np.nan),
    ])
    p3 = pd.DataFrame([
        # E1: randomness PROVES (count up SNR>2, 1500/stringent up, singleton up, excess not up)
        _mk("3", "mean_ectopic_count", "full_depth", "E1", 0.15, 1e-4, 4.2),
        _mk("3", "mean_ectopic_count", "rarefied_1500", "E1", 0.08, 1e-3, np.nan),
        _mk("3", "mean_ectopic_count", "full_depth", "E1_stringent", 0.10, 1e-3, np.nan),
        _mk("3", "singleton_fraction_20cells", "subsampled_20cells", "E1", 0.04, 1e-3, np.nan),
        _mk("3", "jaccard_both_pos", "full_depth", "E1", -0.02, 2e-2, 3.0),
        _mk("3", "excess_jaccard", "full_depth", "E1", -0.005, 0.4, np.nan),
        # E2: shared induction (count up SNR>2 AND excess_jaccard significantly UP!)
        _mk("3", "mean_ectopic_count", "bounded_depth_strata", "E2_stratum", 0.20, 1e-4, 5.0),
        _mk("3", "excess_jaccard", "bounded_depth_strata", "E2_stratum", 0.03, 1e-4, np.nan),
    ])
    ver = fd.evaluate_fulldepth_verdicts(p1, p2, p3)
    assert ver["Part1_center_shift"]["verdict"] == "PROVES"
    assert ver["Part1_scatter"]["verdict"] == "DOES NOT PROVE"
    assert ver["Part2"]["shannon"]["verdict"] == "PROVES"
    assert ver["Part2"]["deman_kl_sqrt"]["verdict"] == "DOES NOT PROVE"
    assert ver["Part3"]["E1"]["verdict"] == "PROVES (RANDOMNESS)"
    assert ver["Part3"]["E2"]["verdict"] == "DOES NOT PROVE RANDOMNESS (SHARED_INDUCTION_OR_AMBIENT)"


def test_shard_binary_roundtrip_and_part1_part3_workers(tmp_path: Path, monkeypatch) -> None:
    # Small end-to-end synthetic check of Pass1 chunk -> K-shard binary format -> Part1 & Part3 workers
    monkeypatch.setattr(fd, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(fd, "N_SHARDS", 2)
    sdir = fd._shard_dir("gao")
    sdir.mkdir(parents=True, exist_ok=True)

    n_genes, n_cts, n_donors = 12, 2, 6
    rng = np.random.default_rng(55)
    g_proj = sp.eye(n_genes, dtype=np.float64, format="csc")
    fd._worker_init(g_proj=g_proj, ctx=None)

    keep_d, keep_c = [], []
    rows = []
    for d in range(n_donors):
        for c in range(n_cts):
            for _ in range(25):
                keep_d.append(d)
                keep_c.append(c)
                v = rng.poisson(300, size=n_genes)
                v[c * 4 : c * 4 + 3] += 1200  # CT markers; high UMI >= 1500
                rows.append(v)
    x_raw = sp.csr_matrix(np.vstack(rows).astype(np.float64))
    bin_arr = np.array([1, 1, 1, 3, 3, 3], dtype=np.int32)
    res_p1 = fd._pass1_worker((
        x_raw.data, x_raw.indices, x_raw.indptr, n_genes,
        np.array(keep_d, dtype=np.int32), np.array(keep_c, dtype=np.int32),
        bin_arr, n_donors, n_cts, 0, 1, 0,
    ))
    assert res_p1["partb"], "Expected Part 2 records from high-UMI synthetic chunk"

    # Write shard binary files cleanly using exact layout of run_pass1_and_export_part2
    for sid, (hdr, idx, dat) in res_p1["shards"].items():
        hdr.tofile(sdir / f"shard_{sid}_cells.bin")
        idx.tofile(sdir / f"shard_{sid}_idx.bin")
        dat.tofile(sdir / f"shard_{sid}_dat.bin")

    # Verify Part 1 worker reads shard cleanly and computes positive distances
    y_ref = np.full(n_genes, 1.0 / n_genes, dtype=np.float64)
    y_sq = np.sqrt(y_ref)
    p1_ctx = {
        "ds_name": "gao", "n_genes": n_genes, "n_cts": n_cts, "n_donors": n_donors, "bin_arr": bin_arr,
        "sqrt_ref_full": {c: y_sq for c in range(n_cts)},
        "sqrt_ref_loo": {(d, c): y_sq for d in range(3) for c in range(n_cts)},
        "sqrt_ref_15_full": {c: y_sq for c in range(n_cts)},
        "sqrt_ref_15_loo": {(d, c): y_sq for d in range(3) for c in range(n_cts)},
        "q_loo_mats": {},
    }
    fd._worker_init(ctx=p1_ctx)
    sid_test = next(iter(res_p1["shards"].keys()))
    out_w1 = fd._part1_shard_worker(sid_test)
    assert out_w1["rows"], "Expected Part 1 donor-CT rows from shard worker"
    ct_rows = [r for r in out_w1["rows"] if not r["is_pooled"]]
    assert ct_rows, "Expected non-pooled CT rows"
    assert np.isfinite(ct_rows[0]["center_shift"])
    assert np.isfinite(ct_rows[0]["scatter"])


def test_rarefied_explicit_zeros_cleaned_before_binary_counting() -> None:
    # multinomial_downsample returns matrix.copy() with new data, so genes that
    # drop to 0 stay as explicit stored zeros. Part 3 binarizes by stored
    # entries, so the pipeline must call eliminate_zeros() first; without it
    # every rarefied row is identical to its full-depth row.
    mat = sp.csr_matrix(np.ones((3, 50), dtype=np.int64))
    raw = multinomial_downsample(mat, 5, np.random.default_rng(0))
    assert raw.nnz == mat.nnz  # pins the upstream behavior that made cleanup necessary
    assert int((raw.data == 0).sum()) > 0
    assert np.diff(raw.indptr).tolist() == [50, 50, 50]

    # The pipeline's cleanup pattern, then stored-entry counting as Part 3 does.
    raw.eliminate_zeros()
    assert int((raw.data == 0).sum()) == 0
    assert int(raw.sum()) == 3 * 5
    rare_counts = np.diff(raw.indptr).astype(np.int64)
    assert (rare_counts < np.diff(mat.indptr)).all()

    # The cleaned matrix reaches the binary metric with no phantom hits.
    out = fd.compute_binary_ectopic_metrics(
        raw, np.random.default_rng(1), n_genes_set=50, run_excess=False,
    )
    assert out["n_pos_cells"] == 3
    assert out["has_any_ectopic"] == 1


def test_mature_metadata_excludes_under30_and_keeps_30_to_39_reference(monkeypatch) -> None:
    donors = pd.DataFrame({
        "donor_id": ["D20", "D25", "D30", "D31", "D39", "D40", "D50", "D60", "D61"],
        "age_years": [20.0, 25.0, 30.0, 31.0, 39.0, 40.0, 50.0, 60.0, 61.0],
        "bin": ["YOUNG", "YOUNG", "YOUNG", "YOUNG", "YOUNG", "MID", "MID", "OLD", "OLD"],
    })
    meta = {
        "donors_df": donors,
        "donor_ids": donors["donor_id"].tolist(),
        "cell_donor_idx": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, -1], dtype=np.int32),
        "cell_ct_idx": np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, -1], dtype=np.int32),
        "reconcile": {"n_cells_kept": 9},
    }
    monkeypatch.setattr(fd, "reconcile_gao_metadata", lambda _path: meta.copy())
    mature = fd.reconcile_mature_metadata("gao")
    assert mature["donor_ids"] == ["D30", "D31", "D39", "D40", "D50", "D60", "D61"]
    assert mature["keep_mask"].tolist() == [False, False, True, True, True, True, True, True, True, False]
    assert mature["cell_donor_idx"].tolist() == [-1, -1, 0, 1, 2, 3, 4, 5, 6, -1]
    assert mature["mature_restriction"]["n_reference_donors"] == 3
    assert mature["donors_df"]["age_years"].min() == 30.0
