#!/usr/bin/env bash
# tools/sonda/S3c_smoke.sh - check smoke di S3c (docs/sonda/REGOLA_S3c.md): come tools/sonda/S3b_smoke.sh (test, catena ridotta, C-3)
# con il metro appaiato tools/sonda/S3c_D.py (G-ID, C-1p, COV, C-2, C-4, C-5) e identita' con results/sonda/S3c/D.json;
# in piu' controlla nel D.json della prova G-ID, CP-perm e il metodo d'intervallo uguale a quello di riferimento.
export S3_TABLE=results/sonda/S3b/area_values.csv S3_H=results/sonda/S3b/h.csv S3_DREF=results/sonda/S3c/D.json S3_PROVV=results/sonda/S3c/misura_provvisoria.json
export S3_DSCRIPT=tools/sonda/S3c_D.py
bash ~/2025.geo_spatialtrans/tools/sonda/S3_smoke.sh "$@"; rc=$?
cd ~/2025.geo_spatialtrans
.venv/bin/python -B -c '
import json, sys
d = json.load(open("/mnt/micron/geo_spatialtrans/S3_sonda/smoke/D.json")); r = json.load(open("results/sonda/S3c/D.json")); ok = True
for k, v in (("GID", d["GID"]["pass"]), ("CP_perm", d["CP_perm"]["pass"]), ("COV method", d["COV"]["method"] == r["COV"]["method"] and not d["COV"]["method_flag"])):
    print(("ok  " if v else "FAIL"), k); ok &= bool(v)
sys.exit(0 if ok else 1)' || rc=1
[ $rc -eq 0 ] && echo "S3c SMOKE PASS" || echo "S3c SMOKE FAIL"
exit $rc
