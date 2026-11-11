
import json, sys
M = "/mnt/micron/geo_spatialtrans"; ok = 0; tot = 0; out = []
F = ["counts_2um","counts_8um","bins","bins_8um","overlap","overlap_8um","origin","origin_8um","background","background_8um","arrived_2um","cells","genes"]
for t in ["A1_r1","A2_r1","A3_r1","A4_f1","A5_r1","A6_r1"]:
    a = json.load(open(f"{M}/S2.8/final504/{t}_k0.md5")); b = json.load(open(f"{M}/S2.7/p1/new/{t}_k0.md5"))
    pairs = [(f, a[f], b[f]) for f in F] + [("layers(S2.7 names)", a["layers_s27"], b["layers"])]
    for f, x, y in pairs:
        tot += 1; ok += x == y; out.append(f"{t},{f},{x},{y},{x == y}")
open("/home/user/2025.geo_spatialtrans/results/S2.8/S2.8_G2bis_final_code.csv", "w").write("tag,field,md5_final,md5_S2.7,equal\n" + "\n".join(out) + "\n")
print(f"G-2 bis: {ok}/{tot} impronte uguali a S2.7 con il codice finale")
