#!/usr/bin/env python
"""tools/R3b_consensus.py — R3b: maschera CO «nuclei di Space Ranger con la forma di StarDist» (linfonodo A4, cuore A6).

Pre-registrazione: docs/steps/R3b.md (cancelli 48573f5, ok di Luca ea1e22a). Definizioni operative (dal file dello step):
 - appaiamento reciproco sulla sovrapposizione in pixel fra etichette Space Ranger (SR) e StarDist (SD): (i, j) appaiati se
   j = argmax_j ov(i, j), i = argmax_i ov(i, j), ov > 0; pareggi all'identificativo piu' piccolo;
 - raster: SR senza corrispondente dipinti per primi, SD appaiati sopra (vincono i pixel SD); SR rimasti senza pixel tolti;
   etichette nuove 1..n in ordine di etichetta SR; i nuclei visti solo da StarDist non entrano;
 - tabella dei nuclei: tools/R2_common.nuclei_table con le maschere di tessuto e bolle di R2 (flag keep ricalcolati sulla CO);
 - filtro cp2 (P-2): etichette con centroide nel tessuto valido, come tools/R2c_candidates.cp2_lab.
Mutanti (variabile R3B_MUTANT): M1 SD solo-StarDist inclusi; M2 SR dipinti sopra SD; M3 appaiamento non reciproco (miglior SD
per ogni SR); M4 SR rimasti vuoti non tolti; M5 filtro cp2 saltato.
Uso (root del repo): .venv/bin/python tools/R3b_consensus.py build [ncores]   -> $R3B_OUT/{masks,nuclei,regress}
                     .venv/bin/python tools/R3b_consensus.py check [ncores]   -> $R3B_OUT/check/G4_<A>_<roi>.csv (cancello G-4)
"""
import os, sys
import numpy as np, pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
MUT = os.environ.get("R3B_MUTANT", "none")
OUT = Path(os.environ.get("R3B_OUT", "/mnt/micron/geo_spatialtrans/R3b"))
MASKS = Path("/mnt/micron/geo_spatialtrans/R2/masks")
TISSUES = ("A4", "A6")


def overlap_table(sr, sd):
    """Sovrapposizioni in pixel fra etichette > 0 di due maschere della stessa forma: DataFrame (sr, sd, ov)."""
    k = (sr > 0) & (sd > 0)
    a = sr[k].astype(np.int64); b = sd[k].astype(np.int64)
    if a.size == 0:
        return pd.DataFrame({"sr": np.zeros(0, np.int64), "sd": np.zeros(0, np.int64), "ov": np.zeros(0, np.int64)})
    M = int(b.max()) + 1
    u, n = np.unique(a * M + b, return_counts=True)
    return pd.DataFrame({"sr": u // M, "sd": u % M, "ov": n.astype(np.int64)})


def mutual_match(ov):
    """Coppie (sr, sd) di massima sovrapposizione reciproca; pareggi all'identificativo piu' piccolo."""
    if len(ov) == 0:
        return ov.copy()
    best_sd = ov.sort_values(["sr", "ov", "sd"], ascending=[True, False, True]).drop_duplicates("sr")
    if MUT == "M3":
        return best_sd.sort_values("sr").reset_index(drop=True)
    best_sr = ov.sort_values(["sd", "ov", "sr"], ascending=[True, False, True]).drop_duplicates("sd")
    m = best_sd.merge(best_sr[["sr", "sd"]], on=["sr", "sd"])
    return m.sort_values("sr").reset_index(drop=True)


def build_consensus(sr, sd):
    """Maschera CO da due maschere di etichette (stessa forma). Restituisce (etichette int64, tabella per etichetta nuova)."""
    sr = np.asarray(sr, np.int64); sd = np.asarray(sd, np.int64)
    assert sr.shape == sd.shape
    ov = overlap_table(sr, sd); m = mutual_match(ov)
    sr_ids = np.unique(sr[sr > 0]); sd_ids = np.unique(sd[sd > 0])
    sr_only = np.setdiff1d(sr_ids, m["sr"].to_numpy())
    nsr = int(sr.max()) + 1; nsd = int(sd.max()) + 1
    lut_sr = np.zeros(nsr, np.int64); lut_sr[sr_only] = sr_only              # etichetta provvisoria = etichetta SR
    lut_sd = np.zeros(nsd, np.int64); lut_sd[m["sd"].to_numpy()] = m["sr"].to_numpy()
    sd_only = np.setdiff1d(sd_ids, m["sd"].to_numpy())
    if MUT == "M1":
        lut_sd[sd_only] = nsr + np.arange(len(sd_only))
    pa = lut_sr[sr]; pb = lut_sd[sd]
    if MUT == "M2":
        prov = np.where(pa > 0, pa, pb)
    else:
        prov = np.where(pb > 0, pb, pa)
    present = np.unique(prov[prov > 0])
    cand = np.union1d(np.union1d(m["sr"].to_numpy(), sr_only), (nsr + np.arange(len(sd_only))) if MUT == "M1" else [])
    keep = cand if MUT == "M4" else present
    keep = np.asarray(keep, np.int64)
    lut_new = np.zeros(int(max(prov.max(), keep.max() if keep.size else 0)) + 1, np.int64)
    lut_new[keep] = np.arange(1, len(keep) + 1)
    lab = lut_new[prov]
    # tabella
    msr = dict(zip(m["sr"].to_numpy(), m["sd"].to_numpy())); mov = dict(zip(m["sr"].to_numpy(), m["ov"].to_numpy()))
    n_sr = np.bincount(sr.ravel(), minlength=nsr); n_sd = np.bincount(sd.ravel(), minlength=nsd)
    n_out = np.bincount(lab.ravel(), minlength=len(keep) + 1)
    rows = []
    for k, p in enumerate(keep, start=1):
        if p >= nsr:
            j = int(sd_only[p - nsr]); rows.append((k, "sd_only", 0, j, 0, 0, int(n_sd[j]), int(n_out[k]))); continue
        j = int(msr.get(p, 0))
        rows.append((k, "sd" if j else "sr", int(p), j, int(mov.get(p, 0)), int(n_sr[p]), int(n_sd[j]) if j else 0, int(n_out[k])))
    tab = pd.DataFrame(rows, columns=["label", "src", "sr_label", "sd_label", "ov_px", "n_px_sr", "n_px_sd", "n_px"])
    tab.attrs.update(n_sr=len(sr_ids), n_sd=len(sd_ids), n_matched=len(m), n_sr_only=len(sr_only),
                     n_sr_only_dropped=int(len(sr_only) - np.isin(sr_only, present).sum()), n_sd_only=len(sd_only))
    return lab, tab


def matched_pixels_identical(co, sd, tab):
    """Cancello G-4: ogni nucleo CO di origine SD coincide pixel per pixel con il suo StarDist, e nessun altro pixel CO cade
    dove c'e' uno StarDist appaiato. Restituisce il numero di pixel discordanti."""
    t = tab[tab.src == "sd"]
    lut = np.zeros(int(sd.max()) + 1, np.int64); lut[t.sd_label.to_numpy()] = t.label.to_numpy()
    expect = lut[np.asarray(sd, np.int64)]
    is_sd = np.zeros(int(co.max()) + 1, bool); is_sd[t.label.to_numpy()] = True
    got = np.where(is_sd[co], co, 0)
    return int(np.sum(expect != got) + np.sum((expect > 0) & (co != expect)))


def cp2_filter(lab, upp, valid, side_um):
    """Etichette con centroide (dei pixel) nel tessuto valido; le altre a 0. Come R2c_candidates.cp2_lab, su un array dato."""
    lab = np.asarray(lab, np.int64).copy()
    if MUT == "M5":
        return lab
    pxv = side_um / valid.shape[1]
    r, c = np.nonzero(lab); l = lab[r, c]; M = int(lab.max()) + 1
    n = np.bincount(l, minlength=M); ids = np.nonzero(n)[0]
    cx = np.bincount(l, c.astype(float), M)[ids] / n[ids] * upp
    cy = np.bincount(l, r.astype(float), M)[ids] / n[ids] * upp
    ci = np.clip((cx / pxv).astype(int), 0, valid.shape[1] - 1); ri = np.clip((cy / pxv).astype(int), 0, valid.shape[0] - 1)
    keep = np.zeros(M, bool); keep[ids[valid[ri, ci]]] = True
    lab[~keep[lab]] = 0
    return lab


def load_native(A, roi, seg):
    return np.load(MASKS / f"{A}_{roi}_{seg}_native.npz")["labels"]


def roi_list():
    import R2_common as C
    r = C.rois_table()
    return [(a, b) for a, b in zip(r.archetype, r.roi_id) if a in TISSUES]


def build_one(A, roi):
    import R2_common as C
    fm = OUT / "masks" / f"{A}_{roi}_consensus_native.npz"
    if fm.exists():
        return "skip"
    sr = load_native(A, roi, "spaceranger"); sd = load_native(A, roi, "stardist_he")
    lab, tab = build_consensus(sr, sd)
    rw = C.rois_table(); rw = rw[(rw.archetype == A) & (rw.roi_id == roi)].iloc[0]; upp = float(rw.um_per_px)
    rgb = C.load_roi(A, roi); t, b = C.masks_for(rgb, upp, 1, A, roi)
    nuc = C.nuclei_table(lab.astype(np.uint32), upp, t, b, A, roi, "consensus", 1)
    nsd = C.nuclei_table(sd, upp, t, b, A, roi, "stardist_he", 1)            # regressione G-3: tabella StarDist di R2
    ref = pd.read_parquet(ROOT / f"results/R2/nuclei/{A}_{roi}_stardist_he_native.parquet")
    cols = [c for c in ref.columns]
    eq = {c: bool(np.array_equal(nsd[c].to_numpy(), ref[c].to_numpy())) for c in cols}
    for d in ("masks", "nuclei", "regress"):
        (OUT / d).mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(archetype=A, roi=roi, n_rows=len(nsd), n_rows_ref=len(ref), n_cols_equal=sum(eq.values()),
                       n_cols=len(cols), cols_diff=";".join(c for c, v in eq.items() if not v))]).to_csv(
        OUT / "regress" / f"nuclei_stardist_{A}_{roi}.csv", index=False)
    tab.to_parquet(OUT / "masks" / f"{A}_{roi}_consensus_map.parquet", index=False)
    pd.DataFrame([dict(archetype=A, roi=roi, **tab.attrs)]).to_csv(OUT / "masks" / f"{A}_{roi}_consensus_counts.csv", index=False)
    nuc.to_parquet(OUT / "nuclei" / f"{A}_{roi}_consensus_native.parquet", index=False)
    np.savez_compressed(str(fm) + ".tmp.npz", labels=lab.astype(np.uint32)); os.replace(str(fm) + ".tmp.npz", fm)
    return "ok"


def _r2c():
    import importlib.util
    spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py")
    K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
    return K


def check_one(A, roi):
    """Cancello G-4 sulla maschera CO salvata, ricalcolato dai file; piu' la regressione del filtro cp2 (G-3)."""
    import hashlib
    K = _r2c(); rois = K.S.rois_table()
    co = np.load(OUT / "masks" / f"{A}_{roi}_consensus_native.npz")["labels"].astype(np.int64)
    tab = pd.read_parquet(OUT / "masks" / f"{A}_{roi}_consensus_map.parquet")
    cnt = pd.read_csv(OUT / "masks" / f"{A}_{roi}_consensus_counts.csv").iloc[0]
    sr = load_native(A, roi, "spaceranger").astype(np.int64); sd = load_native(A, roi, "stardist_he").astype(np.int64)
    labs = np.unique(co[co > 0]); n_px = np.bincount(co.ravel(), minlength=int(co.max()) + 1)
    t_sd = tab[tab.src == "sd"]; t_sr = tab[tab.src == "sr"]
    # pixel degli StarDist appaiati
    lut_sd = np.zeros(int(sd.max()) + 1, np.int64); lut_sd[t_sd.sd_label.to_numpy()] = t_sd.label.to_numpy()
    exp_sd = lut_sd[sd]
    # SR solo-Space Ranger: i suoi pixel CO = pixel SR non coperti da uno StarDist appaiato
    free = (sr > 0) & (exp_sd == 0)
    n_free = np.bincount(sr[free], minlength=int(sr.max()) + 1)
    lut_sr = np.zeros(int(co.max()) + 1, np.int64); lut_sr[t_sr.label.to_numpy()] = t_sr.sr_label.to_numpy()
    k = lut_sr[co] > 0
    lab_c = K.cp2_lab(A, roi, "spaceranger", rois)[0]
    rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]; upp = float(rw.um_per_px); side_um = float(rw.side_px) * upp
    valid = K.S.read_valid(A, roi)
    cp2_sr_eq = bool(np.array_equal(cp2_filter(sr, upp, valid, side_um), lab_c))
    cp2_sd_eq = bool(np.array_equal(cp2_filter(sd, upp, valid, side_um), K.cp2_lab(A, roi, "stardist_he", rois)[0]))
    re_lab, re_tab = build_consensus(sr, sd)
    md5 = lambda a: hashlib.md5(np.ascontiguousarray(a.astype(np.uint32)).tobytes()).hexdigest()
    row = dict(archetype=A, roi=roi, n_labels=len(labs), n_tab=len(tab),
               labels_consecutive=bool(np.array_equal(labs, np.arange(1, len(tab) + 1)) and np.array_equal(labs, tab.label.to_numpy())),
               n_empty=int(np.sum(n_px[labs] == 0)), src_values=";".join(sorted(tab.src.unique())),
               count_identity=bool(len(labs) == int(cnt.n_matched) + int(cnt.n_sr_only) - int(cnt.n_sr_only_dropped)),
               sd_pixels_bad=matched_pixels_identical(co, sd, tab), sd_pixels_bad_indep=int(np.sum((exp_sd > 0) != (np.isin(co, t_sd.label.to_numpy())))
                                                                                          + np.sum((exp_sd > 0) & (co != exp_sd))),
               sr_pixels_bad=int(np.sum(sr[k] != lut_sr[co][k])),
               sr_npx_bad=int(np.sum(n_free[t_sr.sr_label.to_numpy()] != n_px[t_sr.label.to_numpy()])),
               cp2_sr_equal=cp2_sr_eq, cp2_sd_equal=cp2_sd_eq,
               md5_saved=md5(co), md5_rerun=md5(re_lab), rerun_identical=bool(np.array_equal(co, re_lab)))
    (OUT / "check").mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(OUT / "check" / f"G4_{A}_{roi}.csv", index=False)
    return "ok"


def _run(job):
    stage, A, roi = job
    try:
        return f"[{stage}] {A} {roi} " + (build_one(A, roi) if stage == "build" else check_one(A, roi))
    except Exception as e:  # noqa: BLE001
        return f"[{stage}] {A} {roi} ERRORE: {e!r}"


if __name__ == "__main__":
    from multiprocessing import Pool
    stage = sys.argv[1]; nc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    jobs = [(stage, a, r) for a, r in roi_list()]
    with Pool(nc) as p:
        res = p.map(_run, jobs, chunksize=1)
    print("\n".join(res)); bad = [x for x in res if "ERRORE" in x]
    print(f"[{stage}] {len(res)} unita', {len(bad)} errori"); sys.exit(1 if bad else 0)
