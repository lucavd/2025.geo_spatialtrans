"""tools/R5_lumens.py — R5, CP-2 (BL-086): lumi interni vuoti all'H&E, UMI osservati contro previsti dal nucleo dei bordi esterni.
Uso: .venv/bin/python tools/R5_lumens.py <KEY> <outdir> <params.json>   (params: K2 mediano della sezione, scarto, a, b dal file delle stime)
Tessere di 150 x 150 bin con in_tissue di Space Ranger >= 0.5, a >= 100 µm dal bordo dell'area, dentro l'immagine; regione = tessera + 75 bin.
Lumo = componente di vetro H&E sulla sottogriglia, >= 2500 µm², che non tocca il bordo della regione, bolle <= 1 % dell'area, centroide nella tessera.
Bin del lumo: interamente di vetro; fasce di distanza dal bordo del lumo 2-4, 4-8, 8-16, 16-32 µm. Scrive <outdir>/lumens/<KEY>_lumens.csv."""
import sys, os, json, time
import numpy as np, pandas as pd
from scipy import ndimage
sys.path.insert(0, "tools")
import R5_common as C
import R5_tiles as T5

BANDS = [(2, 4), (4, 8), (8, 16), (16, 32)]

if __name__ == "__main__":
    KEY, OUT, PJ = sys.argv[1], sys.argv[2], sys.argv[3]; ds = T5.DS[KEY]
    prm = json.load(open(PJ))[KEY]
    os.makedirs(f"{OUT}/lumens", exist_ok=True)
    t0 = time.time()
    z = np.load(f"{OUT}/section/{KEY}.npz"); U, IT, M, t = z["U"], z["in_tissue"], z["M"], z["t"]
    nr, nc = U.shape
    p, H, W, upp = T5.image_info(ds)
    PXR, PXC = T5.pixel_grids(ds, nr, nc)
    TB, MB = T5.TB, T5.MB
    rows_out = []
    n_tiles = 0
    for i in range(nr // TB):
        for j in range(nc // TB):
            R0, C0 = i * TB, j * TB
            if IT[R0:R0 + TB, C0:C0 + TB].mean() < 0.5: continue
            if R0 < T5.BORDER_BINS or C0 < T5.BORDER_BINS or R0 + TB > nr - T5.BORDER_BINS or C0 + TB > nc - T5.BORDER_BINS: continue
            cr = [(M[0, 0] * r + M[0, 1] * c + t[0], M[1, 0] * r + M[1, 1] * c + t[1]) for r in (R0 - 0.5, R0 + TB - 0.5) for c in (C0 - 0.5, C0 + TB - 0.5)]
            if not all(0 <= a < H and 0 <= b < W for a, b in cr): continue
            n_tiles += 1
            Rr, Cr, nb = R0 - MB, C0 - MB, TB + 2 * MB
            Ml, tl_, _ = T5.local_affine(PXR, PXC, Rr, Cr, nb)
            f, fa, bub, ok = T5.region_masks(p, H, W, upp, Ml, tl_, Rr, Cr, nb)
            glass = (~f) & ok
            lab, nl = ndimage.label(glass)
            if nl == 0: continue
            area = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, nl + 1)) * C.SUB_UM ** 2
            border = set(np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])) - {0}
            bubf = ndimage.mean(bub, lab, index=np.arange(1, nl + 1))
            com = ndimage.center_of_mass(glass, lab, index=np.arange(1, nl + 1))
            cand = [k + 1 for k in range(nl) if area[k] >= 2500 and (k + 1) not in border and bubf[k] <= 0.01
                    and MB * 4 <= com[k][0] < (MB + TB) * 4 and MB * 4 <= com[k][1] < (MB + TB) * 4]
            if not cand: continue
            Ureg = np.zeros((nb, nb)); inarr = np.zeros((nb, nb), bool)
            a0, a1, b0, b1 = max(Rr, 0), min(Rr + nb, nr), max(Cr, 0), min(Cr + nb, nc)
            Ureg[a0 - Rr:a1 - Rr, b0 - Cr:b1 - Cr] = U[a0:a1, b0:b1]; inarr[a0 - Rr:a1 - Rr, b0 - Cr:b1 - Cr] = True
            rows, cols = np.indices(Ureg.shape); rows += Rr; cols += Cr
            P = C.prepare_region(Ureg, f.astype(float), rows, cols)
            fit = {"sigma_r": prm["sigma"], "sigma_c": prm["sigma"], "w": prm["w"], "ell": prm["ell"], "dr": prm["dr"], "dc": prm["dc"], "a": prm["a"], "b": prm["b"]}
            lam = C.predict(fit, P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], P["c"])
            inimg = C.bin_occupancy(ok.astype(float)) >= 1
            for k in cand:
                lm = (lab == k).astype(float)
                occ_l = C.bin_occupancy(lm)                       # frazione del bin nel lumo
                dl = ndimage.distance_transform_edt(lab == k) * C.SUB_UM
                d_bin = dl[C.SUB // 2::C.SUB, C.SUB // 2::C.SUB][:nb, :nb]
                full = (occ_l >= 1) & inarr & inimg
                rec = {"key": KEY, "tile": f"t{R0}_{C0}", "lumen": int(k), "area_um2": float(area[k - 1]), "n_bins_full": int(full.sum())}
                for lo, hi in BANDS:
                    m = full & (d_bin >= lo) & (d_bin < hi)
                    rec[f"n_{lo}_{hi}"] = int(m.sum()); rec[f"obs_{lo}_{hi}"] = float(Ureg[m].sum()); rec[f"pred_{lo}_{hi}"] = float(lam[m].sum())
                rows_out.append(rec)
    L = pd.DataFrame(rows_out); L.to_csv(f"{OUT}/lumens/{KEY}_lumens.csv", index=False)
    print(f"[{KEY}] interior tiles {n_tiles} lumens {len(L)} {time.time()-t0:.0f}s", flush=True)
