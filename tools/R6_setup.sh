#!/usr/bin/env bash
# R6 setup: python deps, marker atlas, output dirs
set -euo pipefail
cd ~/2025.geo_spatialtrans
mkdir -p /mnt/micron/geo_spatialtrans/R6 results/R6
~/.local/bin/uv pip install --python .venv/bin/python scanpy leidenalg igraph scikit-image shapely 2>&1 | tail -3
.venv/bin/python -c "import scanpy, leidenalg, skimage, shapely; print('scanpy', scanpy.__version__, 'leidenalg', leidenalg.__version__)"
M=/mnt/micron/geo_spatialtrans/R6/PanglaoDB_markers_27_Mar_2020.tsv.gz
[ -s "$M" ] || curl -sSL -o "$M" https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz
ls -la "$M"; md5sum "$M"; zcat "$M" | head -2
