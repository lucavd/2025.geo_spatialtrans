"""tools/R4_test.py — cancello G-1 di R4: batteria di test su dati sintetici e 3 mutanti.
Uso (root del repo): .venv/bin/python tools/R4_test.py   -> PASS/FAIL per asserzione; exit 1 se un FAIL."""
import sys, os, importlib
import numpy as np, scipy.sparse as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import R4_common as C

RES = []
def check(name, cond):
    RES.append((name, bool(cond))); return bool(cond)

def suite(M, verbose=True):
    """Esegue la batteria contro il modulo M; ritorna la lista (nome, esito)."""
    out = []
    def ck(n, c): out.append((n, bool(c)))
    # 1. aggregazione 2 -> 8: griglia 8 x 8 con valori noti
    rr, cc = np.indices((8, 8)); v = (rr * 8 + cc).ravel()
    agg = M.aggregate_totals(rr.ravel(), cc.ravel(), v)
    exp = {(i, j): int(sum(v[(rr.ravel() // 4 == i) & (cc.ravel() // 4 == j)])) for i in range(2) for j in range(2)}
    ck("agg_2to8_values", len(agg) == 4 and all(int(r.umi) == exp[(r.row, r.col)] for r in agg.itertuples()))
    ck("agg_2to8_total_conserved", int(agg.umi.sum()) == int(v.sum()))
    agg2 = M.aggregate_totals(np.array([3, 4]), np.array([0, 0]), np.array([1, 1]))
    ck("agg_2to8_boundary_rows_3_4_split", len(agg2) == 2)
    # 2. strati su un'immagine di etichette sintetica (1 px = 0.25 um): nucleo 1 tenuto, nucleo 2 non tenuto
    upp = 0.25; mask = np.zeros((400, 400), np.int32)
    mask[90:110, 90:110] = 1; mask[300:320, 300:320] = 2
    rows = np.array([100, 100, 100 + int(6 / upp), 100 + int(12 / upp) + 10, 310, 310 + int(4 / upp)])
    cols = np.array([100, 112, 100, 100, 310, 310])
    ln = M.footprint_labels(mask, rows, cols, upp, 8.0, mode="nuc"); le = M.footprint_labels(mask, rows, cols, upp, 8.0, mode="exp")
    s = M.strata(ln, le, np.array([1]))
    ck("strata_codes", list(s) == [1, 2, 2, 0, 3, 3])
    ck("strata_shift_moves_mask", M.shifted_mask(mask, upp, 20.0)[90 + 80, 90 + 80] == 1)
    # 3. stimatore del rapporto nucleo/corona su Poisson simulati
    rng = np.random.default_rng(1)
    st = np.r_[np.ones(200000, int), 2 * np.ones(400000, int)]
    u2 = np.r_[rng.poisson(2.0, 200000), rng.poisson(1.0, 400000)]
    u1 = rng.poisson(1.5, 600000)
    ck("ratio_recovers_2", abs(M.density_ratio(u2, st) - 2.0) < 0.03)
    ck("ratio_recovers_1", abs(M.density_ratio(u1, st) - 1.0) < 0.02)
    ck("ratio_nan_without_ring", np.isnan(M.density_ratio(u1[:10], np.ones(10, int))))
    # 4. P(0) della binomiale negativa contro simulazione
    mu, th = 3.0, 1.5; x = rng.negative_binomial(th, th / (th + mu), 400000)
    ck("nb_zero_formula", abs(M.nb_zero_expected([x.mean()], [x.var()])[0] - (x == 0).mean()) < 0.004)
    ck("nb_zero_poisson_branch", abs(M.nb_zero_expected([2.0], [1.5])[0] - np.exp(-2.0)) < 1e-12)
    # 5. zeri in eccesso: ~0 su binomiale negativa, > 0.02 con 30 % di dropout
    ng, nc = 300, 4000
    means = np.exp(rng.uniform(np.log(2), np.log(30), ng)); thetas = rng.uniform(0.5, 5, ng)
    Xd = np.column_stack([rng.negative_binomial(thetas[g], thetas[g] / (thetas[g] + means[g]), nc) for g in range(ng)])
    X = sp.csr_matrix(Xd)
    genes = M.eligible_genes(X)
    ex = M.excess_zeros(X, genes)
    ck("excess_zero_nb_near_0", len(genes) > 250 and abs(np.median(ex)) < 0.005)
    Xdo = M.apply_dropout(X, 0.30, np.random.default_rng(2))
    ck("dropout_fraction_30pct", abs(1 - Xdo.nnz / X.nnz - 0.30) < 0.005)
    ck("excess_zero_detects_dropout", np.median(M.excess_zeros(Xdo, M.eligible_genes(Xdo))) > 0.02)
    w = np.random.default_rng(3).multinomial(nc, np.ones(nc) / nc)
    m_w, _, _ = M.weighted_moments(X, w); Xr = sp.csr_matrix(np.repeat(Xd, w, axis=0)); m_r, _, _ = M.weighted_moments(Xr)
    ck("weighted_moments_equal_replication", np.allclose(m_w, m_r))
    # 6. bootstrap riproducibile e verdetti
    vals = np.array([1.0, 1.2, 0.9, 1.1, 1.05])
    a = M.boot_one_stage(vals, rng=np.random.default_rng(7)); b = M.boot_one_stage(vals, rng=np.random.default_rng(7))
    ck("boot_reproducible", a[:3] == b[:3])
    g = [np.arange(10.0), np.arange(5.0) + 3]
    t1 = M.boot_two_stage(g, np.median, B=200, rng=np.random.default_rng(7)); t2 = M.boot_two_stage(g, np.median, B=200, rng=np.random.default_rng(7))
    ck("boot2_reproducible", t1[:3] == t2[:3] and t1[1] <= t1[0] <= t1[2])
    ck("verdicts", [M.verdict(0.95, 1.05, 0.909, 1.1), M.verdict(1.2, 1.4, 0.909, 1.1), M.verdict(1.05, 1.2, 0.909, 1.1),
                    M.verdict(1.3, 1.5, 1.25, np.inf), M.verdict(np.nan, 1, 0, 1)] ==
       ["confermata", "smentita", "indeterminata", "confermata", "non valutabile"])
    # 7. spettro: rumore senza picco; rumore con reticolo diagonale -> picco alla frequenza giusta
    img = rng.poisson(3.0, (512, 512)).astype(float); valid = np.ones_like(img, bool)
    S, k = M.welch_spectrum(img, valid); p0 = M.spectral_peak(S)
    yy, xx = np.indices(img.shape); lat = img + 2.0 * (1 + np.cos(2 * np.pi * (yy + xx) / 8.0))
    S1, _ = M.welch_spectrum(lat, valid); p1 = M.spectral_peak(S1)
    ck("spectrum_noise_no_peak", k == 64 and not p0["present"])
    ck("spectrum_lattice_peak", p1["present"] and abs(abs(p1["fy"]) - 0.125) < 0.02 and abs(abs(p1["fx"]) - 0.125) < 0.02)
    # 8. (revisione F01, F10, F16) Nyquist e reticolo diagonale fuori dagli assi
    par = rng.poisson(3.0, (512, 512)).astype(float) * np.where(xx % 2 == 1, 1.5, 1.0)
    Sp, _ = M.welch_spectrum(par, valid); pp = M.spectral_peak(Sp); po = M.spectral_peak(Sp, off_axis=True)
    ck("spectrum_parity_at_nyquist", pp["present"] and abs(abs(pp["fx"]) - 0.5) < 1e-9 and abs(pp["fy"]) < 1e-9)
    ck("spectrum_offaxis_ignores_parity", not po["present"])
    pd_ = par + 2.0 * (1 + np.cos(2 * np.pi * (yy + xx) / 8.0))
    Sd, _ = M.welch_spectrum(pd_, valid); pdia = M.spectral_peak(Sd, off_axis=True)
    ck("spectrum_offaxis_finds_diagonal", pdia["present"] and abs(abs(pdia["fy"]) - 0.125) < 0.02 and abs(abs(pdia["fx"]) - 0.125) < 0.02)
    flat = np.ones((64, 64)); f5 = flat.copy(); f5[32 + 10, 32 + 10] = 5.0; f12 = flat.copy(); f12[32 + 10, 32 + 10] = 12.0
    ck("spectrum_threshold_is_10", (not M.spectral_peak(f5)["present"]) and M.spectral_peak(f12)["present"])
    # 9. rapporto come rapporto di medie (non di mediane)
    ck("ratio_is_ratio_of_means", abs(M.density_ratio(np.array([0, 0, 0, 10, 1, 1, 1, 1]), np.array([1, 1, 1, 1, 2, 2, 2, 2])) - 2.5) < 1e-12)
    # 10. bootstrap: ricampiona davvero le ROI e le cellule; intervallo al 95 %
    v = np.random.default_rng(11).normal(0, 1, 400)
    e1, l1, h1, _ = M.boot_one_stage(v, rng=np.random.default_rng(12))
    se = v.std(ddof=1) / np.sqrt(len(v))
    ck("boot1_interval_95pct", abs((h1 - l1) / (2 * 1.96 * se) - 1) < 0.12)
    gb = [np.full(50, float(k)) for k in range(5)]                 # varianza solo fra ROI
    gw = [np.random.default_rng(20 + k).normal(0, 1, 50) for k in range(5)]   # varianza solo dentro le ROI (ROI quasi uguali)
    gw = [x - x.mean() for x in gw]
    _, lb, hb, _ = M.boot_two_stage(gb, np.mean, B=300, rng=np.random.default_rng(13))
    _, lw, hw, _ = M.boot_two_stage(gw, np.mean, B=300, rng=np.random.default_rng(14))
    ck("boot2_resamples_rois", hb - lb > 0.5)
    ck("boot2_resamples_cells", hw - lw > 0.05)
    # 11. ammissione dei geni per media >= 2 (soglia inclusa)
    Xe = sp.csr_matrix(np.column_stack([np.full(100, 2), np.r_[np.full(99, 2), 1], np.r_[np.full(50, 4), np.zeros(50)]]))
    ck("eligible_genes_by_mean_inclusive", list(M.eligible_genes(Xe)) == [0, 2])
    ck("verdict_boundaries_inclusive", M.verdict(0.909, 1.1, 0.909, 1.1) == "confermata" and M.verdict(1.1, 1.3, 0.909, 1.1) == "indeterminata")
    # 12. (revisione F06) taglie diverse fra cellule: la binomiale negativa semplice vede zeri in difetto, quella con fattore di taglia no
    sfac = np.exp(rng.normal(0, 0.6, nc)); sfac /= sfac.mean()
    Xs = sp.csr_matrix(np.column_stack([rng.negative_binomial(thetas[g], thetas[g] / (thetas[g] + means[g] * sfac), nc) for g in range(ng)]))
    gs = M.eligible_genes(Xs)
    ck("plain_nb_biased_negative_with_size_factors", np.median(M.excess_zeros(Xs, gs)) < -0.005)
    ck("sizefactor_nb_unbiased", abs(np.median(M.nb_sf_excess(Xs[:, gs], sfac))) < 0.003)
    ck("sizefactor_nb_detects_dropout", np.median(M.nb_sf_excess(M.apply_dropout(Xs, 0.3, np.random.default_rng(5))[:, gs], sfac)) > 0.02)
    return out

def mutant(name, attr, fn):
    M = importlib.reload(C); setattr(M, attr, fn)
    res = suite(M); n_fail = sum(1 for _, ok in res if not ok)
    importlib.reload(C)
    return n_fail

if __name__ == "__main__":
    res = suite(C)
    for n, ok in res:
        print(f"{'PASS' if ok else 'FAIL'}  {n}")
    import pandas as pd
    _orig_strata = C.strata
    def m_strata(ln, le, keep):          # strati nucleo e corona scambiati
        s = _orig_strata(ln, le, keep); t = s.copy(); t[s == 1] = 2; t[s == 2] = 1; return t
    def m_agg(rows, cols, umi, f=4):     # indice 2 -> 8 sfasato di uno
        r, c = (np.asarray(rows) + 1) // f, np.asarray(cols) // f
        return pd.DataFrame({"row": r, "col": c, "umi": np.asarray(umi, np.int64)}).groupby(["row", "col"], as_index=False)["umi"].sum()
    def m_nb(mean, var):                  # zeri attesi di Poisson al posto della binomiale negativa
        return np.exp(-np.asarray(mean, float))
    _o1, _o2, _oe, _osp, _over = C.boot_one_stage, C.boot_two_stage, C.eligible_genes, C.spectral_peak, C.verdict
    def m_ratio_median(umi, s):           # rapporto di mediane
        umi = np.asarray(umi, float); s = np.asarray(s); return float(np.median(umi[s == 1]) / max(np.median(umi[s == 2]), 1e-9))
    def m_boot1_noresample(values, stat=np.mean, B=C.B, rng=None):
        e = float(stat(np.asarray(values, float))); return e, e, e, np.full(B, e)
    def m_boot1_90(values, stat=np.mean, B=C.B, rng=None):
        v = np.asarray(values, float); rng = rng if rng is not None else np.random.default_rng(C.SEED)
        reps = np.array([stat(v[i]) for i in rng.integers(0, len(v), size=(B, len(v)))])
        return float(stat(v)), float(np.percentile(reps, 5)), float(np.percentile(reps, 95)), reps
    def m_boot2_noroi(groups, stat, B=C.B, rng=None):
        rng = rng if rng is not None else np.random.default_rng(C.SEED)
        reps = np.array([stat(np.concatenate([g[rng.integers(0, len(g), len(g))] for g in groups])) for _ in range(B)])
        return float(stat(np.concatenate(groups))), float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5)), reps
    def m_boot2_nocell(groups, stat, B=C.B, rng=None):
        rng = rng if rng is not None else np.random.default_rng(C.SEED)
        reps = np.array([stat(np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])) for _ in range(B)])
        return float(stat(np.concatenate(groups))), float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5)), reps
    def m_elig_detect(X, min_mean=C.MIN_MEAN):   # ammissione per frazione di cellule rilevate
        _, _, z = C.weighted_moments(X); return np.flatnonzero(1 - z >= 0.5)
    def m_elig_strict(X, min_mean=C.MIN_MEAN):
        m1, _, _ = C.weighted_moments(X); return np.flatnonzero(m1 > min_mean)
    def m_spec_nonyq(S, r_min=3, thr=10.0, off_axis=False):   # vecchio rilevatore senza Nyquist (F01)
        n = S.shape[0]; cy = n // 2; yy, xx = np.indices(S.shape); rad = np.rint(np.hypot(yy - cy, xx - cy)).astype(int); ratio = np.zeros_like(S)
        for r in range(r_min, cy):
            m = rad == r; md = np.median(S[m]); ratio[m] = S[m] / md if md > 0 else 0
        i = np.unravel_index(np.argmax(ratio), S.shape)
        return {"ratio": float(ratio[i]), "present": bool(ratio[i] > thr), "fy": float((i[0] - cy) / n), "fx": float((i[1] - cy) / n)}
    def m_spec_thr3(S, r_min=3, thr=10.0, off_axis=False): return _osp(S, r_min, 3.0, off_axis)
    def m_verdict_strict(lo, hi, a, b):
        if not (np.isfinite(lo) and np.isfinite(hi)): return "non valutabile"
        if lo > a and hi < b: return "confermata"
        if hi <= a or lo >= b: return "smentita"
        return "indeterminata"
    muts = [("strata_swapped", "strata", m_strata), ("agg_offset", "aggregate_totals", m_agg), ("nb_as_poisson", "nb_zero_expected", m_nb),
            ("ratio_of_medians", "density_ratio", m_ratio_median), ("boot1_no_resampling", "boot_one_stage", m_boot1_noresample),
            ("boot1_90pct", "boot_one_stage", m_boot1_90), ("boot2_no_roi_stage", "boot_two_stage", m_boot2_noroi),
            ("boot2_no_cell_stage", "boot_two_stage", m_boot2_nocell), ("eligible_by_detection", "eligible_genes", m_elig_detect),
            ("eligible_strict", "eligible_genes", m_elig_strict), ("spectrum_without_nyquist", "spectral_peak", m_spec_nonyq),
            ("spectrum_threshold_3", "spectral_peak", m_spec_thr3), ("verdict_strict", "verdict", m_verdict_strict)]
    for nm, attr, fn in muts:
        nf = mutant(nm, attr, fn)
        res.append((f"mutant_{nm}_detected", nf > 0)); print(f"{'PASS' if nf > 0 else 'FAIL'}  mutant_{nm}_detected ({nf} asserzioni fallite)")
    n_fail = sum(1 for _, ok in res if not ok)
    print(f"R4_test: {len(res) - n_fail}/{len(res)} PASS")
    sys.exit(1 if n_fail else 0)
