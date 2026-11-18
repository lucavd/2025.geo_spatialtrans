"""tools/F3.1_PH4d_spill.py - F3.1, seguito di PH-4 nato dopo averlo letto (descrittivo, non pre-registrato): marcatori estranei nel SOLO nucleo
(foreign_nuc, senza allargamento) e dopo l'allargamento (foreign_exp = M-3), vero contro simulato (k = 0, 1) e Q6 (k = 0), con i marcatori estranei dei
profili generati (foreign_profile) e la purezza vera del nucleo (truth_purity_nuc) del simulato. Legge solo i b2c_label.json gia' scritti.
Uso: python tools/F3.1_PH4d_spill.py <main_dir> <out_dir>"""
import sys, glob, json, os, pandas as pd
main, out = sys.argv[1:3]; rows = []
for side, pat in [("real", f"{main}/real/A*"), ("sim", f"{main}/sim/A*_k*"), ("simq6", f"{main}/simq6/A*_k0")]:
    for d in sorted(glob.glob(pat)):
        f = f"{d}/b2c_label.json"
        if not os.path.exists(f): continue
        j = json.load(open(f)); rows.append({"side": side, "area": os.path.basename(d), "A": j["archetype"], **{k: j.get(k) for k in
            ["foreign_nuc", "foreign_exp", "foreign_profile", "truth_purity_nuc", "truth_purity_exp", "truth_foreign_type_nuc", "n_typed_nuc"]}})
R = pd.DataFrame(rows); R.to_csv(f"{out}/F3.1_PH4d_spill_areas.csv", index=False)
S = R.groupby(["A", "side"])[["foreign_nuc", "foreign_exp", "foreign_profile", "truth_purity_nuc", "truth_foreign_type_nuc"]].mean().round(3)
S["n"] = R.groupby(["A", "side"]).size(); S.to_csv(f"{out}/F3.1_PH4d_spill.csv"); pd.set_option("display.width", 200); print(S.to_string())
