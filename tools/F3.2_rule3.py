#!/usr/bin/env python3
"""tools/F3.2_rule3.py - F3.2, terza regola per le aree di verifica (decisione tecnica dell'agente, DEC-095), scritta prima di eseguirla.
Corregge i due difetti della seconda (file dello step): (1) distanza sui punteggi dei marcatori (media su quadratini di 8 um di log1p(conteggi per 10 000)
mediata sui geni dell'insieme), non sui rapporti, che esplodevano per gli insiemi quasi assenti; (2) tolleranza dalla variabilita' fra aree di R2 della
classe: sigma_ks = max(DS fra aree di R2 della classe k, 0.25 * DS fra finestre dello stesso lato su tutta la sezione di riferimento).
Classe k se D_k = max_s |x_s - mu_ks| / sigma_ks <= 3 e D_k < D_j per ogni altra classe del tessuto. Filtri: copertura e grigio della seconda regola;
bianco H&E (tessere di 50 um con piu' di meta' dei pixel a densita' ottica < 0.15) <= massimo delle aree di R2 del tessuto + 0.05 (rene e fegato: <= 0.10).
E' il terzo sguardo sulle stesse 30 aree di R2: la prova non e' cieca. Per questo la scelta finale richiede anche l'accordo di due letture a vista (DEC-095).
Uso: .venv/bin/python tools/F3.2_rule3.py validate            -> results/F3.2/F3.2_rule3_validation.csv
     .venv/bin/python tools/F3.2_rule3.py classify <label>    -> results/F3.2/F3.2_rule3_pool_<label>.csv (classe delle candidate gia' estratte)
     .venv/bin/python tools/F3.2_rule3.py pool500 A4          -> linfonodo con aree di 500 um: candidate, classe, foglio (results/F3.2/sheets/pool500_A4.png)
"""
import sys, math, importlib.util
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
_s = importlib.util.spec_from_file_location("f32areas", str(ROOT / "tools/F3.2_areas.py")); A = importlib.util.module_from_spec(_s); _s.loader.exec_module(A)
OUT = ROOT / "results/F3.2"; Z_MAX = 3.0; FLOOR = 0.25; WHITE_MARGIN = 0.05; WHITE_NONE = 0.10

def white_fraction(img, upp, box):
    t = int(round(A.TILE_UM / upp)); c0, r0 = int(round(box[0] / upp)), int(round(box[1] / upp)); n = int(round((box[2] - box[0]) / A.TILE_UM))
    reg = A.read_region(img, r0, r0 + n * t, c0, c0 + n * t)
    if reg.ndim == 3 and reg.shape[0] in (3, 4) and reg.shape[-1] not in (3, 4): reg = np.moveaxis(reg, 0, -1)
    od = -np.log10((reg[..., :3].astype(np.float32) + 1.0) / 256.0).sum(axis=-1); g = od < 0.15
    H, W = (g.shape[0] // t) * t, (g.shape[1] // t) * t
    return float((g[:H, :W].reshape(H // t, t, W // t, t).mean(axis=(1, 3)) > 0.5).mean())

def xbox(B, box, sets):
    m = (B.x >= box[0]) & (B.x < box[2]) & (B.y >= box[1]) & (B.y < box[3])
    return {s: float(B.loc[m, s].mean()) for s in sets}, float(m.sum() / (((box[2] - box[0]) / 8.0) * ((box[3] - box[1]) / 8.0)))

def section_sd(B, sets, side):
    gx0, gy0 = math.floor(B.x.min() / 50) * 50, math.floor(B.y.min() / 50) * 50
    nx, ny = int((B.x.max() - gx0) // 50) + 1, int((B.y.max() - gy0) // 50) + 1
    ix = ((B.x.values - gx0) // 50).astype(int); iy = ((B.y.values - gy0) // 50).astype(int); n = int(round(side / 50))
    def bs(a):
        S = np.zeros((a.shape[0] + 1, a.shape[1] + 1)); S[1:, 1:] = a.cumsum(0).cumsum(1); return S[n:, n:] - S[:-n, n:] - S[n:, :-n] + S[:-n, :-n]
    cnt = np.zeros((ny, nx)); np.add.at(cnt, (iy, ix), 1); C = bs(cnt); ok = C / ((side / 8.0) ** 2) >= 0.5
    out = {}
    for s in sets:
        a = np.zeros((ny, nx)); np.add.at(a, (iy, ix), B[s].values); out[s] = float((bs(a) / np.maximum(C, 1))[ok].std())
    return out

def refs_r2(tissue, side=1000.0):
    """aree di R2 (o i loro quarti, se side = 500) della sezione di riferimento: punteggi, copertura, bianco"""
    sec = A.REF_SECTION[tissue]; B, u = A.bin_scores(sec); P = A.paths(sec); sets = ["|".join(s) for s in A.set_list(tissue)]
    rows = []
    for a in A.r2_table():
        if a["section"] != sec: continue
        x0, y0, x1, y1 = a["box"]
        boxes = [a["box"]] if side >= 1000 else [(x0 + i * 500, y0 + j * 500, x0 + (i + 1) * 500, y0 + (j + 1) * 500) for j in (0, 1) for i in (0, 1)]
        for q, b in enumerate(boxes):
            x, cov = xbox(B, b, sets)
            rows.append(dict(id=a["id"] + (f"_q{q}" if side < 1000 else ""), parent=a["id"], klass=a["klass"], box=b, x=x, cov=cov, white=white_fraction(P["img"], u, b)))
    return rows, sets, section_sd(B, sets, side)

def classify(x, refs, sets, S, exclude_parent=None):
    D = {}
    for k in sorted({r["klass"] for r in refs}):
        rr = [r for r in refs if r["klass"] == k and r["parent"] != exclude_parent]
        d = 0.0
        for s in sets:
            v = np.array([r["x"][s] for r in rr]); mu = v.mean(); sd = v.std(ddof=1) if len(v) > 1 else 0.0
            d = max(d, abs(x[s] - mu) / max(sd, FLOOR * S[s], 1e-9))
        D[k] = d
    k = min(D, key=D.get); ok = D[k] <= Z_MAX and all(D[k] < D[j] for j in D if j != k)
    return (k if ok else "altro"), D

def validate():
    res = []
    for tissue in ["intestine", "crc", "lymph", "brain", "heart"]:
        refs, sets, S = refs_r2(tissue); cmin = A.cov_min(tissue); P = A.paths(A.REF_SECTION[tissue]); u = A.upp_of(P)
        for r in refs:
            cls, D = classify(r["x"], refs, sets, S, exclude_parent=r["parent"])
            wmax = max(o["white"] for o in refs if o["parent"] != r["parent"]) + WHITE_MARGIN
            grey = A.grey_fraction(P["img"], u, r["box"])
            rec = dict(kind="R2_ROI", id=r["id"], tissue=tissue, class_expected=r["klass"], class_assigned=cls, own_ok=cls == r["klass"],
                       coverage=r["cov"], cov_min=cmin, white=r["white"], white_max=wmax, grey=grey,
                       admitted=(r["cov"] >= cmin and grey <= A.GREY_MAX and r["white"] <= wmax))
            others = [k for k in D if k != r["klass"]]; rec["negatives"] = len(others); rec["negatives_rejected"] = sum(cls != k for k in others)
            for k, v in D.items(): rec[f"D_{k}"] = v
            res.append(rec); print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
    # ritaglio del linfonodo di R1
    P = A.paths("A4"); u = A.upp_of(P); B, _ = A.bin_scores("A4"); refs, _, _ = refs_r2("lymph")
    pos = pd.read_parquet(P["sq2"] / "spatial/tissue_positions.parquet"); it = pos[pos.in_tissue == 1]
    rng = np.random.default_rng(42); s = it.sample(1, random_state=int(rng.integers(1e9))).iloc[0]
    half = int(round(400.0 / u / 2)); r0c, c0c = max(int(s.pxl_row_in_fullres) - half, 0), max(int(s.pxl_col_in_fullres) - half, 0)
    box = (c0c * u, r0c * u, c0c * u + 400.0, r0c * u + 400.0); _, cov = xbox(B, box, []); w = white_fraction(P["img"], u, box)
    wmax = max(o["white"] for o in refs) + WHITE_MARGIN
    res.append(dict(kind="R1_crop", id="A4_crop1", tissue="lymph", class_expected="scartato", coverage=cov, cov_min=A.cov_min("lymph"), white=w, white_max=wmax,
                    grey=A.grey_fraction(P["img"], u, box), admitted=(cov >= A.cov_min("lymph") and w <= wmax), negatives=0, negatives_rejected=0))
    d = pd.DataFrame(res); d.to_csv(OUT / "F3.2_rule3_validation.csv", index=False); r = d[d.kind == "R2_ROI"]
    print(f"terza regola: (i) classe propria {int(r.own_ok.sum())}/30; (ii) negativi {int(r.negatives_rejected.sum())}/{int(r.negatives.sum())}; "
          f"(iii) ammesse {int(r.admitted.sum())}/30; (iv) ritaglio di R1 scartato {not bool(d[d.kind=='R1_crop'].admitted.iloc[0])}")
    print(r.groupby("tissue").own_ok.sum().to_dict())

def classify_pool(label, side=1000.0, pool_file=None):
    P = A.paths(label); tissue = P["tissue"]; B, u = A.bin_scores(label)
    pf = pool_file or OUT / f"F3.2_pool_{label}.csv"; pool = pd.read_csv(pf); pool = pool[pool.level != "test_tile"]
    if tissue == "none":
        refs, sets, S, wmax = [], [], {}, WHITE_NONE
    else:
        refs, sets, S = refs_r2(tissue, side); wmax = max(r["white"] for r in refs) + WHITE_MARGIN
    rows = []
    for _, c in pool.iterrows():
        box = (c.x0, c.y0, c.x1, c.y1); w = white_fraction(P["img"], u, box)
        if tissue == "none": cls, D = "tessuto", {}
        else:
            x, _ = xbox(B, box, sets); cls, D = classify(x, refs, sets, S)
        rec = dict(label=label, order=int(c.order), area_id=c.area_id, x0=c.x0, y0=c.y0, x1=c.x1, y1=c.y1, side_um=side, rule3_class=cls if w <= wmax else "scartata_bianco",
                   white=w, white_max=wmax, coverage=c.coverage, grey=c.grey_frac)
        for k, v in D.items(): rec[f"D_{k}"] = v
        rows.append(rec)
    o = pd.DataFrame(rows); fo = OUT / (f"F3.2_rule3_pool_{label}.csv" if side >= 1000 else f"F3.2_rule3_pool500_{label}.csv"); o.to_csv(fo, index=False)
    print(label, o.rule3_class.value_counts().to_dict(), "->", fo, flush=True)

def pool500(label="A4"):
    """linfonodo: candidate di 500 um (stessa estrazione della selezione a vista, lato 500, generatore [seme, sezione, 5]); classe con i quarti di R2"""
    A.SIDE_UM = 500.0
    P = A.paths(label); tissue = P["tissue"]; B, u = A.bin_scores(label); img = P["img"]
    refs, sets, S = refs_r2(tissue, 500.0); cmin = min(r["cov"] for r in refs) - 0.05; wmax = max(r["white"] for r in refs) + WHITE_MARGIN
    used = pd.read_csv(OUT / "F3.2_used_areas.csv")
    blocked = [tuple(v) for v in used[(used.section == label) & used.set.isin(["R2_ROI", "R2b_WINDOW"])][["x0", "y0", "x1", "y1"]].values]
    cx, cy = float(np.median(B.x)), float(np.median(B.y)); blocked.append((cx - 125, cy - 125, cx + 125, cy + 125))
    gx0, gy0 = math.floor(B.x.min() / 50) * 50, math.floor(B.y.min() / 50) * 50
    nx, ny = int((B.x.max() - gx0) // 50) + 1, int((B.y.max() - gy0) // 50) + 1
    ix = ((B.x.values - gx0) // 50).astype(int); iy = ((B.y.values - gy0) // 50).astype(int); n = 10
    cnt = np.zeros((ny, nx)); np.add.at(cnt, (iy, ix), 1)
    Sx = np.zeros((ny + 1, nx + 1)); Sx[1:, 1:] = cnt.cumsum(0).cumsum(1); C = (Sx[n:, n:] - Sx[:-n, n:] - Sx[n:, :-n] + Sx[:-n, :-n]) / ((500 / 8.0) ** 2)
    cand = np.argwhere(C >= cmin); rng = np.random.default_rng([A.SEED, sorted(list(A.CAL) + list(A.HOLD)).index(label), 5]); cand = cand[rng.permutation(len(cand))]
    rows, taken = [], []
    for j, i in cand:
        box = (gx0 + i * 50, gy0 + j * 50, gx0 + i * 50 + 500, gy0 + j * 50 + 500)
        if any(A.gap(box, b) < A.MARGIN_UM for b in blocked + taken): continue
        if A.grey_fraction(img, u, box) > A.GREY_MAX or white_fraction(img, u, box) > wmax: continue
        taken.append(box); rows.append(dict(label=label, level="L1", order=len(taken), area_id=f"{label}_h{len(taken):02d}", x0=box[0], y0=box[1], x1=box[2], y1=box[3],
                                            coverage=float(C[j, i]), cov_min=cmin, grey_frac=0.0, um_per_px=u, n_admissible=int(len(cand))))
        if len(taken) == 20: break
    pf = OUT / f"F3.2_pool500_{label}.csv"; pd.DataFrame(rows).to_csv(pf, index=False)
    print(label, "500 um: ammessi", len(cand), "candidate", len(taken), flush=True)
    classify_pool(label, 500.0, pf)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    refs_img = [r for r in refs]; d = pd.read_csv(pf); ncol = 5; nrow = 2 + math.ceil(len(d) / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.25 * nrow), facecolor="white"); axs = np.atleast_2d(axs)
    for ax in axs.ravel(): ax.set_axis_off()
    def thumb(box):
        t = int(round(500 / u)); c0, r0 = int(round(box[0] / u)), int(round(box[1] / u)); reg = A.read_region(img, r0, r0 + t, c0, c0 + t)
        if reg.ndim == 3 and reg.shape[0] in (3, 4) and reg.shape[-1] not in (3, 4): reg = np.moveaxis(reg, 0, -1)
        s = t // 2 * 2; return reg[:s, :s, :3].reshape(s // 2, 2, s // 2, 2, 3).mean(axis=(1, 3)).astype(np.uint8)
    for ri, k in enumerate(["follicolo", "paracorticale"]):
        for ci, r in enumerate([r for r in refs_img if r["klass"] == k][:5]):
            axs[ri, ci].imshow(thumb(r["box"])); axs[ri, ci].set_title(f"RIFERIMENTO {k.upper()}\n{r['id']}", fontsize=8, color="darkgreen")
    for k, (_, r) in enumerate(d.iterrows()):
        ax = axs[2 + k // ncol, k % ncol]; ax.imshow(thumb((r.x0, r.y0, r.x1, r.y1))); ax.set_title(f"#{k + 1}", fontsize=12, fontweight="bold")
    fig.suptitle("A4 (aree di 500 um): per ogni candidata, a quale riga somiglia? (follicolo / paracorticale / nessuna)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98)); fig.savefig(OUT / "sheets" / "pool500_A4.png", dpi=80); plt.close(fig)

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "validate": validate()
    elif cmd == "classify": classify_pool(sys.argv[2])
    elif cmd == "pool500": pool500(sys.argv[2] if len(sys.argv) > 2 else "A4")
