"""tools/sonda/S8_to_misura.py - S8 (docs/sonda/REGOLA_S8.md, «Uscite»): copia A5 e il suo errore standard da results/sonda/S8/A5.json in
results/misura.json (affermazioni.A5), senza ricalcolo. Lo chiama check smoke dopo S3c_smoke.sh."""
import json, os
R = os.path.expanduser("~/2025.geo_spatialtrans/results")
a = json.load(open(f"{R}/sonda/S8/A5.json")); m = json.load(open(f"{R}/misura.json"))
m.setdefault("affermazioni", {})["A5"] = {"valore": a["A5"], "se": a["A5_se"]}
json.dump(m, open(f"{R}/misura.json", "w"), indent=1); print("ok   results/misura.json A5 =", round(a["A5"], 4), "se", round(a["A5_se"], 4), "(da results/sonda/S8/A5.json)")
