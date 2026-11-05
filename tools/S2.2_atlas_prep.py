"""S2.2: prepara le cellule degli atlanti a sonde (Chromium Flex) per la gara: matrice cellule x geni, tipo del catalogo del preset,
sottotipo. Uso: .venv/bin/python tools/S2.2_atlas_prep.py crc | ln | brain
crc: Oliveira 2025 (doi:10.1038/s41588-025-02193-3), paziente P2CRC, QCFilter == Keep; tipo da Level1 (tabella pre-registrata),
sottotipo = Level2 dentro il tipo. ln, brain: senza annotazione; tipi con la pipeline di R6 (Leiden + criterio corretto PH-1).
Scrive /mnt/micron/geo_spatialtrans/S2.2/atlas/<nome>_cells.parquet, _X.npz (cellule x geni), _genes.tsv.
Pre-registrazione: docs/steps/S2.2.md (sezione «Il banco della gara»)."""
import sys, os, numpy as np, pandas as pd, scipy.sparse as sp, h5py
A = "/mnt/micron/geo_spatialtrans/S2.2/atlas"
NAME = sys.argv[1]
L1_MAP = {"Tumor": "epiteliale", "Intestinal Epithelial": "epiteliale", "Fibroblast": "fibroblasti", "Endothelial": "endoteliale",
          "T cells": "immunitario", "B cells": "immunitario", "Myeloid": "immunitario"}   # Smooth Muscle, Neuronal: fuori dai preset


def read_h5(path):
    with h5py.File(path, "r") as h:
        m = h["matrix"]; shape = tuple(int(v) for v in m["shape"][:])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape)   # geni x cellule
        bc = m["barcodes"][:].astype(str); ids = m["features"]["id"][:].astype(str); names = m["features"]["name"][:].astype(str)
        ft = m["features"]["feature_type"][:].astype(str)
    g = ft == "Gene Expression"
    return X[g].T.tocsr(), bc, ids[g], names[g]


if NAME == "crc":
    X, bc, ids, names = read_h5(f"{A}/crc/crc_count_filtered_feature_bc_matrix.h5")
    meta = pd.read_csv(f"{A}/crc/SingleCell_MetaData.csv.gz").set_index("Barcode").reindex(bc)
    sel = ((meta.Patient == "P2CRC") & (meta.QCFilter == "Keep")).to_numpy()
    meta = meta[sel].copy(); X = X[sel]
    meta["cell_type"] = meta.Level1.map(L1_MAP)
    keep = meta.cell_type.notna().to_numpy()
    meta = meta[keep]; X = X[keep]
    cells = pd.DataFrame({"barcode": meta.index, "cell_type": meta.cell_type.to_numpy(), "subtype_raw": meta.Level2.astype(str).to_numpy(),
                          "umi": np.asarray(X.sum(axis=1)).ravel()})
else:
    raise SystemExit(f"atlante {NAME}: link non ancora disponibili")
sp.save_npz(f"{A}/{NAME}_X.npz", X.astype(np.int32))
cells.to_parquet(f"{A}/{NAME}_cells.parquet", index=False)
pd.DataFrame({"gene_id": ids, "gene_name": names}).to_csv(f"{A}/{NAME}_genes.tsv", sep="\t", index=False)
print(f"[{NAME}] cells {X.shape[0]} genes {X.shape[1]}; umi median {np.median(cells.umi):.0f}")
print(cells.cell_type.value_counts().to_string())
