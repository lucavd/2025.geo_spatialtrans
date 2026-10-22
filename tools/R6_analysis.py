"""R6 analysis: per-type references, predictions P-1..P-6, counter-proofs CP-1, CP-2, sensitivities.
Usage: python tools/R6_analysis.py   (reads /mnt/micron/geo_spatialtrans/R6/labels, writes results/R6/*.csv)"""
import os, re, sys, json, numpy as np, pandas as pd
VARIANT = sys.argv[1] if len(sys.argv) > 1 else "prereg"   # prereg | ph1
SUF = "" if VARIANT == "prereg" else "_" + VARIANT
M = "/mnt/micron/geo_spatialtrans/R6"; OUT = "results/R6"; os.makedirs(OUT, exist_ok=True)
SEED = 20261022; B = 2000; rng = np.random.default_rng(SEED)
PRIMARY = ("exp", 20, 20261022)
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
refs = pd.read_csv("results/S1.5/S1.5_real_refs_roi.csv")
refs["area_net_mm2"] = refs.n / refs.density_seg_mm2
NONTYPE = {"non classificabile", "non risolto"}

def load_labels(A, fp="exp", umi=20, seed=20261022):
    f = f"{M}/labels/{A}_{fp}_umi{umi}_seed{seed}_nuclei{SUF}.parquet"
    if not os.path.exists(f): return None
    lab = pd.read_parquet(f)
    # territories (S1.4/real): join on (x, y) of the generator = nucleus centroid
    parts = []
    for roi in lab.roi_id.unique():
        t = pd.read_parquet(f"{M}/terr/{A}_{roi}.parquet"); t["roi_id"] = roi; parts.append(t)
    t = pd.concat(parts, ignore_index=True)
    t["kx"] = t.x.round(4); t["ky"] = t.y.round(4); lab["kx"] = lab.x_um.round(4); lab["ky"] = lab.y_um.round(4)
    lab = lab.merge(t[["roi_id", "kx", "ky", "a_t", "a_n", "interior", "matched", "ecc_n"]], on=["roi_id", "kx", "ky"], how="left")
    lab["nt"] = np.where(lab.interior & lab.matched & (lab.a_n > 0), lab.a_n / lab.a_t, np.nan)
    return lab

def per_roi_type(lab):
    rows = []
    for (roi, ty), g in lab.groupby(["roi_id", "cell_type"]):
        n_roi = (lab.roi_id == roi).sum(); A = lab.archetype.iloc[0]
        area = refs.query("archetype == @A and roi_id == @roi").area_net_mm2.iloc[0]
        rows.append(dict(archetype=A, roi_id=roi, cell_type=ty, n=len(g), n_roi=n_roi, fraction=len(g) / n_roi, density_mm2=len(g) / area,
                         median_area_nuc=g.area_um2.median(), median_ecc_nuc=g.eccentricity.median(),
                         median_a_t_interior=g.a_t[g.interior == True].median(), median_nt=g.nt.median(), n_interior=int((g.interior == True).sum()),
                         median_umi=g.umi.median()))
    return pd.DataFrame(rows)

def boot_stat(lab, fn, B=B, two_stage=True, rois=None):
    """block bootstrap over ROIs (unit), optionally resampling nuclei within each resampled ROI; fn(df) -> float."""
    groups = {r: g for r, g in lab.groupby("roi_id")}; keys = list(groups) if rois is None else [r for r in rois if r in groups]
    est = fn(pd.concat([groups[k] for k in keys])); vals = np.empty(B)
    for b in range(B):
        pick = rng.choice(keys, size=len(keys), replace=True)
        parts = [groups[k].sample(frac=1.0, replace=True, random_state=int(rng.integers(1 << 31))) if two_stage else groups[k] for k in pick]
        vals[b] = fn(pd.concat(parts))
    vals = vals[np.isfinite(vals)]
    return est, np.quantile(vals, 0.025) if len(vals) else np.nan, np.quantile(vals, 0.975) if len(vals) else np.nan, len(vals)

def verdict(lo, hi, lower=None, upper=None):
    """interval of relevance [lower, upper] (None = unbounded); confirmed if CI fully inside, refuted if fully outside, else indeterminate."""
    lower = -np.inf if lower is None else lower; upper = np.inf if upper is None else upper
    if not (np.isfinite(lo) and np.isfinite(hi)): return "non valutabile"
    if lo >= lower and hi <= upper: return "confermata"
    if hi < lower or lo > upper: return "smentita"
    return "indeterminata"

GLIA = {"astrociti", "oligodendrociti", "microglia"}
def ratio_neuron_glia(df):
    a = df.area_um2[df.cell_type == "neuroni"].median(); g = df.area_um2[df.cell_type.isin(GLIA)].median(); return a / g
def ratio_max_min(df, min_frac=0.05):
    d = df[~df.cell_type.isin(NONTYPE)]; fr = d.cell_type.value_counts(normalize=True); ty = fr[fr >= min_frac].index
    med = d[d.cell_type.isin(ty)].groupby("cell_type").area_um2.median(); return med.max() / med.min() if len(med) >= 2 else np.nan
def ecc_diff_cm(df):
    d = df[~df.cell_type.isin(NONTYPE)]; return d.eccentricity[d.cell_type == "cardiomiociti"].median() - d.eccentricity[d.cell_type != "cardiomiociti"].median()
def frac_of(ty):
    return lambda df: (df.cell_type == ty).mean()

# ---- Space Ranger mapping (P-6)
def sr_map(A, coarse, fine):
    s = f"{coarse} | {fine}".lower()
    if A == "A6":
        if "cardiac muscle" in s or "cardiocyte" in s or "cardiac myocyte" in s or "myoblast" in s or coarse.lower() == "muscle cell": return "cardiomiociti"
        if "smooth muscle" in s or "pericyte" in s or "mural" in s: return "muscolo liscio e periciti"
        if "endothelial" in s: return "endoteliali"
        if "fibroblast" in s or "stromal" in s or "connective" in s: return "fibroblasti"
        if any(k in s for k in ["macrophage", "b cell", "t cell", "lymphocyte", "mast", "monocyte", "granulocyte", "plasma cell", "lymphoid", "hematopoietic", "leukocyte", "dendritic"]): return "immunitari"
        return "altro"
    if A == "A5":
        if "astrocyte" in s: return "astrociti"
        if "oligodendrocyte" in s: return "oligodendrociti"
        if "microglia" in s: return "microglia"
        if "neuron associated" in s: return "altro"
        if "neuron" in s or "interneuron" in s or "pyramidal" in s: return "neuroni"
        if "endothelial" in s or "mural" in s or "pericyte" in s or "smooth muscle" in s: return "vascolari"
        return "altro"
    return "altro"
SRFILE = {"A5": "/mnt/micron/geo_spatialtrans/data_real/Visium_HD_6p5mm_Mouse_Brain/segmented_outputs_x/segmented_outputs/cell_types/10x_Cloud/cell_types.csv",
          "A6": "/mnt/micron/geo_spatialtrans/data_real/Visium_HD_6p5mm_Human_Heart/segmented_outputs_x/segmented_outputs/cell_types/10x_Cloud/cell_types.csv"}
def add_sr(lab, A):
    if A not in SRFILE: return lab
    ct = pd.read_csv(SRFILE[A]); ct["sr_type"] = [sr_map(A, c, f) for c, f in zip(ct.coarse_cell_type.astype(str), ct.fine_cell_type.astype(str))]
    ct["sr_lowumi"] = ct.coarse_cell_type.astype(str).str.contains("Low UMI") | (ct.coarse_cell_type.astype(str) == "cell")
    lab = lab.merge(ct[["barcode", "coarse_cell_type", "fine_cell_type", "sr_type", "sr_lowumi"]], left_on="sr_cell_id", right_on="barcode", how="left")
    return lab

def agreement_all(df):
    d = df[df.sr_cell_id.notna() & df.sr_type.notna() & ~df.sr_lowumi.fillna(True)]; return (d.cell_type == d.sr_type).mean()
def agreement_informative(df):
    d = df[df.sr_cell_id.notna() & df.sr_type.notna() & ~df.sr_lowumi.fillna(True) & (df.sr_type != "altro") & ~df.cell_type.isin(NONTYPE)]; return (d.cell_type == d.sr_type).mean()

# ================= main
labs = {}; prt = []
for A in ARCH:
    lab = load_labels(A, *PRIMARY); lab = add_sr(lab, A); labs[A] = lab; prt.append(per_roi_type(lab))
    print(A, "joined territories:", lab.a_t.notna().mean().round(4), "n", len(lab), flush=True)
prt = pd.concat(prt, ignore_index=True); prt.to_csv(f"{OUT}/R6_refs_roi_type{SUF}.csv", index=False)

# pooled per archetype x type with ROI bootstrap (fraction, density; two-stage for medians)
rows = []
for A, lab in labs.items():
    for ty in sorted(lab.cell_type.unique()):
        fr = boot_stat(lab, frac_of(ty), two_stage=False)
        dens = prt.query("archetype == @A and cell_type == @ty").density_mm2
        sub = lab[lab.cell_type == ty]
        ar = boot_stat(lab, lambda d, ty=ty: d.area_um2[d.cell_type == ty].median()) if len(sub) >= 20 else (sub.area_um2.median(), np.nan, np.nan, 0)
        ec = boot_stat(lab, lambda d, ty=ty: d.eccentricity[d.cell_type == ty].median()) if len(sub) >= 20 else (sub.eccentricity.median(), np.nan, np.nan, 0)
        at = boot_stat(lab, lambda d, ty=ty: d.a_t[(d.cell_type == ty) & (d.interior == True)].median()) if (sub.interior == True).sum() >= 20 else (np.nan, np.nan, np.nan, 0)
        nt = boot_stat(lab, lambda d, ty=ty: d.nt[d.cell_type == ty].median()) if sub.nt.notna().sum() >= 20 else (np.nan, np.nan, np.nan, 0)
        rows.append(dict(archetype=A, cell_type=ty, n=len(sub), n_rois_present=sub.roi_id.nunique(), fraction=fr[0], fraction_lo=fr[1], fraction_hi=fr[2],
                         density_mm2_mean_roi=dens.mean(), density_mm2_min_roi=dens.min(), density_mm2_max_roi=dens.max(),
                         median_area_nuc=ar[0], area_lo=ar[1], area_hi=ar[2], median_ecc_nuc=ec[0], ecc_lo=ec[1], ecc_hi=ec[2],
                         median_a_t_interior=at[0], a_t_lo=at[1], a_t_hi=at[2], median_nt=nt[0], nt_lo=nt[1], nt_hi=nt[2], median_umi=sub.umi.median()))
    print(A, "types done", flush=True)
pooled = pd.DataFrame(rows); pooled.to_csv(f"{OUT}/R6_refs_type{SUF}.csv", index=False)

# ---- predictions
SPREAD = {"A1": 1.095, "A2": 1.161, "A3": 1.216, "A4": 1.062, "A5": 2.096, "A6": 1.808}
P = []
e = boot_stat(labs["A5"], ratio_neuron_glia); P.append(dict(id="P-1", archetype="A5", quantity="median nuclear area neurons / glia", estimate=e[0], lo=e[1], hi=e[2], lower=2.1, upper=np.nan, verdict=verdict(e[1], e[2], lower=2.1)))
for A in ARCH:
    e = boot_stat(labs[A], ratio_max_min); P.append(dict(id="P-2", archetype=A, quantity="max/min median nuclear area among types >= 5 %", estimate=e[0], lo=e[1], hi=e[2], lower=SPREAD[A], upper=np.nan, verdict=verdict(e[1], e[2], lower=SPREAD[A])))
e = boot_stat(labs["A6"], ecc_diff_cm); P.append(dict(id="P-3", archetype="A6", quantity="median nuclear eccentricity cardiomyocytes - others", estimate=e[0], lo=e[1], hi=e[2], lower=0.104, upper=np.nan, verdict=verdict(e[1], e[2], lower=0.104)))
l4 = labs["A4"]; fol = ["f1", "f2"]; par = ["p1", "p2", "p3"]
def bdiff(df):
    return (df.cell_type[df.roi_id.isin(fol)] == "linfociti B").mean() - (df.cell_type[df.roi_id.isin(par)] == "linfociti B").mean()
# bootstrap: resample ROIs within each compartment separately (two-stage)
gf = {r: g for r, g in l4[l4.roi_id.isin(fol)].groupby("roi_id")}; gp = {r: g for r, g in l4[l4.roi_id.isin(par)].groupby("roi_id")}
vals = []
for b in range(B):
    pf = [gf[k].sample(frac=1, replace=True, random_state=int(rng.integers(1 << 31))) for k in rng.choice(list(gf), len(gf))]
    pp = [gp[k].sample(frac=1, replace=True, random_state=int(rng.integers(1 << 31))) for k in rng.choice(list(gp), len(gp))]
    vals.append(bdiff(pd.concat(pf + pp)))
P.append(dict(id="P-4", archetype="A4", quantity="B fraction follicular (f1,f2) - paracortical (p1-p3)", estimate=bdiff(l4), lo=np.quantile(vals, .025), hi=np.quantile(vals, .975), lower=0.20, upper=np.nan, verdict=verdict(np.quantile(vals, .025), np.quantile(vals, .975), lower=0.20)))
DOM = {"A1": ("epiteliale", 0.88, None), "A2": ("epiteliale tumorale", 0.70, None), "A3": (["fibroblasti", "muscolo liscio"], 0.75, None), "A4": ("linfociti B", 0.80, fol), "A5": ("neuroni", 0.40, None), "A6": ("cardiomiociti", 0.85, None)}
for A, (ty, hyp, rois) in DOM.items():
    fn = (lambda df, ty=ty: df.cell_type.isin(ty).mean()) if isinstance(ty, list) else frac_of(ty)
    e = boot_stat(labs[A], fn, two_stage=False, rois=rois)
    P.append(dict(id="P-5", archetype=A, quantity=f"fraction of {ty if isinstance(ty, str) else '+'.join(ty)}" + (" (follicular ROIs)" if rois else ""), estimate=e[0], lo=e[1], hi=e[2], lower=hyp - 0.15, upper=hyp + 0.15, verdict=verdict(e[1], e[2], lower=hyp - 0.15, upper=hyp + 0.15), design_value=hyp))
for A in ["A5", "A6"]:
    e = boot_stat(labs[A], agreement_informative, two_stage=False); P.append(dict(id="P-6", archetype=A, quantity="agreement with Space Ranger 10x Cloud coarse type (matched, informative)", estimate=e[0], lo=e[1], hi=e[2], lower=0.75, upper=np.nan, verdict=verdict(e[1], e[2], lower=0.75)))
    e = boot_stat(labs[A], agreement_all, two_stage=False); P.append(dict(id="P-6 (all matched)", archetype=A, quantity="agreement incl. SR types outside the catalogue and our unresolved", estimate=e[0], lo=e[1], hi=e[2], lower=0.75, upper=np.nan, verdict=verdict(e[1], e[2], lower=0.75)))
P = pd.DataFrame(P); P.to_csv(f"{OUT}/R6_predictions{SUF}.csv", index=False); print(P.to_string())

# confusion tables vs Space Ranger
for A in ["A5", "A6"]:
    d = labs[A]; d = d[d.sr_cell_id.notna() & ~d.sr_lowumi.fillna(True)]
    pd.crosstab(d.cell_type, d.sr_type).to_csv(f"{OUT}/R6_sr_confusion_{A}{SUF}.csv")
    print(A, "matched to SR:", len(d), "of", len(labs[A]))

# ---- CP-1: permutation null within ROI (200 permutations)
CP = []
for A in ARCH:
    lab = labs[A]; r_null = []; e_null = []
    for b in range(200):
        sh = lab.copy(); sh["cell_type"] = sh.groupby("roi_id").cell_type.transform(lambda s: s.sample(frac=1, random_state=int(rng.integers(1 << 31))).values)
        r_null.append(ratio_max_min(sh))
        if A == "A6": e_null.append(ecc_diff_cm(sh))
        if A == "A5": e_null.append(ratio_neuron_glia(sh))
    CP.append(dict(id="CP-1", archetype=A, quantity="max/min median area, labels permuted within ROI", median_null=np.nanmedian(r_null), q975_null=np.nanquantile(r_null, .975), threshold=1.05, verdict="PASS" if np.nanquantile(r_null, .975) <= 1.05 else "FAIL"))
    if A == "A6": CP.append(dict(id="CP-1", archetype=A, quantity="ecc diff cardiomyocytes - others, permuted", median_null=np.median(e_null), q975_null=np.quantile(np.abs(e_null), .975), threshold=0.01, verdict="PASS" if np.quantile(np.abs(e_null), .975) <= 0.01 else "FAIL"))
    if A == "A5": CP.append(dict(id="CP-1", archetype=A, quantity="neuron/glia area ratio, permuted", median_null=np.nanmedian(e_null), q975_null=np.nanquantile(e_null, .975), threshold=1.05, verdict="PASS" if np.nanquantile(e_null, .975) <= 1.05 else "FAIL"))
# ---- CP-2: shifted footprint
for A in ARCH:
    TC = "type" if VARIANT == "prereg" else "type_ph1"
    cl_p = pd.read_csv(f"{M}/labels/{A}_exp_umi20_seed20261022_clusters{SUF}.csv"); f_sh = f"{M}/labels/{A}_shift_umi20_seed20261022_clusters{SUF}.csv"
    if not os.path.exists(f_sh): continue
    cl_s = pd.read_csv(f_sh)
    row = dict(id="CP-2", archetype=A, quantity="fraction of clusters unresolved: primary vs shifted footprint", primary=(cl_p[TC] == "non risolto").mean(), shifted=(cl_s[TC] == "non risolto").mean())
    if A in SRFILE:
        ls = add_sr(load_labels(A, "shift", 20, 20261022), A); row["agreement_shift"] = agreement_informative(ls); row["agreement_primary"] = agreement_informative(labs[A])
        row["verdict"] = "PASS" if row["agreement_shift"] < 0.5 else "FAIL"
    CP.append(row)
CP = pd.DataFrame(CP); CP.to_csv(f"{OUT}/R6_counterproofs{SUF}.csv", index=False); print(CP.to_string())

# ---- sensitivities: nuc footprint, umi 10/40, second seed -> label agreement and ARI with primary
from sklearn.metrics import adjusted_rand_score
S = []
for A in ARCH:
    base = labs[A].set_index(["roi_id", "label"])
    for fp, umi, seed, name in [("nuc", 20, 20261022, "footprint nucleus only"), ("exp", 10, 20261022, "UMI >= 10"), ("exp", 40, 20261022, "UMI >= 40"), ("exp", 20, 7, "Leiden seed 7")]:
        alt = load_labels(A, fp, umi, seed)
        if alt is None: continue
        alt = alt.set_index(["roi_id", "label"]).reindex(base.index)
        both = base.cell_type.notna() & alt.cell_type.notna() & ~base.cell_type.isin(NONTYPE) & ~alt.cell_type.isin(NONTYPE)
        S.append(dict(archetype=A, variant=name, n_both_typed=int(both.sum()), label_agreement=(base.cell_type[both] == alt.cell_type[both]).mean(),
                      ari_clusters=adjusted_rand_score(base.cluster[both].astype(str), alt.cluster[both].astype(str)), frac_nontype_alt=alt.cell_type.isin(NONTYPE).mean()))
S = pd.DataFrame(S); S.to_csv(f"{OUT}/R6_sensitivity{SUF}.csv", index=False); print(S.to_string())
json.dump({"B": B, "seed": SEED, "primary": PRIMARY, "variant": VARIANT}, open(f"{OUT}/R6_analysis_params{SUF}.json", "w"))
