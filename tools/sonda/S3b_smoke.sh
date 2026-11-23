#!/usr/bin/env bash
# tools/sonda/S3b_smoke.sh - check smoke di S3b (docs/sonda/REGOLA_S3b.md): lo stesso di tools/sonda/S3_smoke.sh (test, catena ridotta, C-3)
# con il metro riparato: tabella e h di results/sonda/S3b/, identita' con results/sonda/S3b/D.json.
export S3_TABLE=results/sonda/S3b/area_values.csv S3_H=results/sonda/S3b/h.csv S3_DREF=results/sonda/S3b/D.json S3_PROVV=results/sonda/S3b/misura_provvisoria.json
exec bash ~/2025.geo_spatialtrans/tools/sonda/S3_smoke.sh "$@"
