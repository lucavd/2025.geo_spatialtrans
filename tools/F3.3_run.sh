#!/usr/bin/env bash
# tools/F3.3_run.sh - F3.3: catena sulle aree di verifica (L1, L2) e sulle aree di taratura del tumore del colon con i geni comuni.
# Fasi: segmentazione (GPU), StarDist e consenso (linfonodo, cuore), regioni, esportazione, bin2cell, RCTD (W processi), misure.
# Rilanciabile: salta cio' che esiste. Uso: bash tools/F3.3_run.sh [W]
set -uo pipefail
cd ~/2025.geo_spatialtrans
W=${1:-4}; M=/mnt/micron/geo_spatialtrans; R=$M/F3.3/runs; REF=$M/F3.1/main/ref; mkdir -p $R $M/F3.3/measures
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
.venv/bin/python tools/F3.3_segment.py holdout > $M/F3.3/log_segment.log 2>&1 || echo "SEGMENT FAIL"
tools/py_stardist.sh tools/F3.3_stardist.py holdout > $M/F3.3/log_stardist.log 2>&1 || echo "STARDIST FAIL"
.venv/bin/python tools/F3.3_consensus.py holdout > $M/F3.3/log_consensus.log 2>&1 || echo "CONSENSUS FAIL"
Rscript --vanilla tools/F3.3_regions.R holdout > $M/F3.3/log_regions.log 2>&1 || echo "REGIONS FAIL"
NONCRC=$(.venv/bin/python -c "import pandas as pd; a=pd.read_csv('results/F3.3/F3.3_areas.csv'); print(','.join(a[(a.level!='IN')&(~a.crc)].area))")
CRCH=$(.venv/bin/python -c "import pandas as pd; a=pd.read_csv('results/F3.3/F3.3_areas.csv'); print(','.join(a[(a.level!='IN')&(a.crc)].area))")
CRCI=$(.venv/bin/python -c "import pandas as pd; a=pd.read_csv('results/F3.3/F3.3_areas.csv'); print(','.join(a[(a.level=='IN')&(a.crc)].area))")
.venv/bin/python tools/F3.3_export.py panel $NONCRC > $M/F3.3/log_export.log 2>&1 || echo "EXPORT FAIL panel"
.venv/bin/python tools/F3.3_export.py common $CRCH,$CRCI >> $M/F3.3/log_export.log 2>&1 || echo "EXPORT FAIL common"
one() {
  local d=$1; local A=$(python3 -c "import json;print(json.load(open('$d/export.json'))['archetype'])")
  [ -f $d/b2c_label.json ] || .venv-b2c/bin/python tools/F3.1_b2c.py $d $REF $A real label > $d/b2c_label.log 2>&1 || { echo "B2C FAIL $d"; return; }
  [ -f $d/rctd.json ] || Rscript --vanilla tools/F3.1_rctd.R $d $REF $A real > $d/rctd.log 2>&1 || { echo "RCTD FAIL $d"; return; }
  local o=$M/F3.3/measures/$(basename $d).json
  [ -s $o ] || .venv-b2c/bin/python tools/F3.3_holdout_measures.py $d $o > $o.log 2>&1 || echo "MEASURE FAIL $d"
  echo "done $d"
}
export -f one; export M REF
ls -d $R/*/ | sed 's#/$##' > $M/F3.3/jobs_runs.txt
xargs -a $M/F3.3/jobs_runs.txt -P $W -L 1 bash -c 'one $0' > $M/F3.3/log_runs.log 2>&1
grep -c '^done' $M/F3.3/log_runs.log; grep -v '^done' $M/F3.3/log_runs.log | head
