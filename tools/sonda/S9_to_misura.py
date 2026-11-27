"""tools/sonda/S9_to_misura.py - S9 (docs/sonda/REGOLA_S9.md, «Esito e conseguenze»): copia in results/misura.json, senza ricalcolo,
A2 (results/sonda/S9/A2_L.json), A4 = e_U (results/sonda/S9/A4_U.json; esito con C in A4_C.json) e A5 col guadagno
(results/sonda/S9/A5_gain.json, decisione di Luca 2026-11-26), che sostituisce il coseno di S8. Lo chiama check smoke dopo S8_to_misura.py."""
import json, os
R = os.path.expanduser("~/2025.geo_spatialtrans/results"); S = f"{R}/sonda/S9"
a2 = json.load(open(f"{S}/A2_L.json")); u = json.load(open(f"{S}/A4_U.json"))["tumour+stroma"]; c = json.load(open(f"{S}/A4_C.json"))["tumour+stroma"]
a5 = json.load(open(f"{S}/A5_gain.json")); m = json.load(open(f"{R}/misura.json")); af = m.setdefault("affermazioni", {})
af["A2"] = {"valore": a2["A2"], "se": a2["A2_se"]}
af["A4"] = {"valore": u["e"], "se": u["e_se"], "catena_C": {"valore": c["e"], "se": c["e_se"], "esito": c["esito"]}}
af["A5"] = {"valore": a5["A5"], "se": a5["A5_se"], "punteggio": "guadagno sulla verita permutata"}
json.dump(m, open(f"{R}/misura.json", "w"), indent=1)
print("ok   results/misura.json A2 =", a2["A2"], "se", round(a2["A2_se"], 4), "| A4 (U) =", round(u["e"], 4), "se", round(u["e_se"], 4), "| A4 (C) =", round(c["e"], 4), "| A5 (guadagno) =", round(a5["A5"], 4), "se", round(a5["A5_se"], 4))
