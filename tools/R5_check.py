"""tools/R5_check.py — G-3: lettura dei dati. (a) grezza (pannello della filtrata) = filtrata sui bin in tessuto; (b) 8 µm = somma dei 16 bin
a 2 µm; (c) totali per cellula e strati = R4 (dalle tappe del rapporto); (d) coordinate dei bin = quelle di R4 (stessi file: md5 uguali a R4).
Uso: .venv/bin/python tools/R5_check.py <outdir>. Scrive results/R5/R5_G3.txt."""
import sys, json, glob
O = sys.argv[1]; lines = []; ok_all = True
for k in ("A1", "A23", "A4", "A5", "A6"):
    d = json.load(open(f"{O}/section/{k}_G3.json"))
    a = d["G3a_raw_equals_filtered_in_tissue"] and d["G3a_filtered_all_in_tissue"]; b = d["G3b_n_8um_bins_differing"] == 0
    lines.append(f"{k} (a) {'PASS' if a else 'FAIL'} [come scritta, tutte le feature della grezza: {'PASS' if d['G3a_as_written_all_raw_features_equal'] else 'FAIL'}, "
                 f"{d['G3a_as_written_n_bins_differing']} bin diversi; grezza {d['n_features_raw']} feature, filtrata {d['n_features_filtered']}] (b) {'PASS' if b else 'FAIL'}")
    ok_all &= a and b
c = {}
for f in glob.glob(f"{O}/nc/A?_*__K2.json"):
    j = json.load(open(f)); c[j["tag"]] = j.get("G3c_identical_to_R4", False)
nc_ok = len(c) == 30 and all(c.values())
lines.append(f"(c) totali per cellula identici a R4 in {sum(c.values())}/{len(c)} ROI: {'PASS' if nc_ok else 'FAIL'}")
r4 = {l.split()[1]: l.split()[0] for l in open("results/R4/R4_md5_before.txt") if "tissue_positions.parquet" in l and "square_002um" in l}
r5 = {l.split()[1]: l.split()[0] for l in open("results/R5/R5_md5_before.txt") if "tissue_positions.parquet" in l and "square_002um" in l}
same = [r4.get(p) == m for p, m in r5.items() if p in r4]
dd = len(same) > 0 and all(same)
lines.append(f"(d) posizioni dei bin a 2 µm con md5 uguali a quelle lette in R4: {sum(same)}/{len(same)} {'PASS' if dd else 'FAIL'}")
ok_all &= nc_ok and dd
lines.append(f"G-3: {'PASS' if ok_all else 'FAIL'}")
open("results/R5/R5_G3.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
