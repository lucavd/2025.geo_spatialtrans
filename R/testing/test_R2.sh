#!/usr/bin/env bash
# R/testing/test_R2.sh — asserzioni della sessione R2 dalle tabelle in results/R2 (un comando, PASS/FAIL per asserzione).
# Non rilancia la segmentazione (vedi reports/R2.qmd, sez. Riproducibilità); verifica coerenza, riproducibilità e soglie pre-registrate.
set -uo pipefail
cd "$(dirname "$0")/../.."
.venv/bin/python - <<'EOF'
import pandas as pd, numpy as np, sys
R = "results/R2/"; fails = 0
def chk(name, ok, detail=""):
    global fails; print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}"); fails += (not ok)
rois = pd.read_csv(R+"R2_rois_checked.csv"); summ = pd.read_csv(R+"R2_roi_summary.csv"); nat = summ[summ.scale==1]
rep = pd.read_csv(R+"R2_repro.csv"); srk = pd.read_csv(R+"R2_spaceranger_check.csv"); cons = pd.read_csv(R+"R2_consensus_density.csv")
nuc = pd.read_parquet(R+"R2_nuclei_all.parquet")
chk("C-R2.4 30 ROI dentro l'immagine", len(rois)==30 and rois.inside_image.all(), f"n={len(rois)}")
chk("C-R2.4 copertura in_tissue > 0.7 in tutti i ROI", (rois.cov_in_tissue > 0.7).all(), f"min={rois.cov_in_tissue.min():.2f}")
chk("C-R2.2 run ripetuti bit-identici", len(rep)==4 and rep.identical_labels.all(), f"{rep.identical_labels.sum()}/{len(rep)}")
k = nuc[nuc.keep & (nuc.scale==1)].head(100000)
chk("C-R2.3 area_um2 = eq_diam consistente", np.allclose(k.eq_diam_um, 2*np.sqrt(k.area_um2/np.pi)), "")
chk("C-R2.3 um/px per dataset in [0.2737, 0.2741]", nuc[nuc.scale==1].um_per_px.between(0.2737, 0.2741).all(), "")
chk("C-R2.5 OD ematossilina dentro > anello in >= 95% dei poligoni SR (25 ROI)", (srk.frac_od_inside_gt_ring >= 0.95).all(), f"min={srk.frac_od_inside_gt_ring.min():.3f}")
chk("C-R2.6 maschera bolle applicata in A6 (frazione riportata > 0)", (nat[(nat.archetype=="A6")&(nat.method=="cellpose_rgb")].bubble_frac > 0).all(), "")
chk("C-R2.7 frazione oggetti flaggati < 5%", ((nat.frac_small + nat.frac_large) < 0.05).all(), f"max={(nat.frac_small+nat.frac_large).max():.3f}")
chk("completezza: 30 ROI x 3 scale x {cellpose_rgb, stardist_he}", len(summ[summ.method.isin(["cellpose_rgb","stardist_he"])])==180, f"n={len(summ[summ.method.isin(['cellpose_rgb','stardist_he'])])}")
chk("completezza: spaceranger su 25 ROI (A2-A6)", len(nat[nat.method=="spaceranger"])==25, "")
sp = pd.read_csv(R+"spatial/R2_spatial_summary.csv")
chk("finestra spatstat coerente: nessun nucleo keep fuori finestra (30 ROI)", (sp.n_lost_window == 0).all() and len(sp)==30, f"max persi={sp.n_lost_window.max()}")
cA = cons.groupby("archetype").density_consensus.median()
print("--- attese B (pre-registrate; un FAIL qui e' un risultato, non un bug):")
chk("B-R2.2 A4 densita' consenso in [1900, 6300]", 1900 <= cA["A4"] <= 6300, f"{cA['A4']:.0f}")
chk("B-R2.4 A5 densita' consenso in [1900, 3500]", 1900 <= cA["A5"] <= 3500, f"{cA['A5']:.0f}")
chk("B-R2.5 A6 densita' consenso in [600, 2000]", 600 <= cA["A6"] <= 2000, f"{cA['A6']:.0f}")
chk("CP-R2.5 A3 e A2 distinguibili (intervalli disgiunti)", cons[cons.archetype=="A3"].density_consensus.max() < cons[cons.archetype=="A2"].density_consensus.min(), "")
print(f"\n{fails} FAIL"); sys.exit(0)
EOF
