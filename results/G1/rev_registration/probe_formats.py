import json, h5py, pyarrow.parquet as pq, numpy as np
D="/mnt/micron/geo_spatialtrans/data_lung/X1"
pf=pq.ParquetFile(f"{D}/xenium/Xenium_V1_Human_Lung_Cancer_FFPE_transcripts.parquet")
print(pf.schema_arrow); print("rows", pf.metadata.num_rows)
print(pq.ParquetFile(f"{D}/xenium/Xenium_V1_Human_Lung_Cancer_FFPE_cells.parquet").schema_arrow.names)
pn=json.load(open(f"{D}/xenium/Xenium_V1_Human_Lung_Cancer_FFPE_gene_panel.json"))
print(pn.keys(), pn["payload"].keys()); t=pn["payload"]["targets"]; print(len(t)); print(json.dumps(t[0])[:600])
B=f"{D}/visium/binned/binned_outputs"
with h5py.File(f"{B}/square_002um/filtered_feature_bc_matrix.h5") as f:
    g=f["matrix"]; print(list(g.keys()), list(g["features"].keys()), g["shape"][:], g["data"].shape, g["data"].dtype)
    print(g["features/id"][:3], g["features/name"][:3], set(x.decode() for x in g["features/feature_type"][:]))
tp=pq.read_table(f"{B}/square_002um/spatial/tissue_positions.parquet").to_pandas(); print(tp.shape); print(tp.head(3).to_string()); print(tp.in_tissue.value_counts().to_dict())
print(tp[["array_row","array_col"]].max().to_dict())
tp8=pq.read_table(f"{B}/square_008um/spatial/tissue_positions.parquet").to_pandas(); print(tp8.shape, tp8[["array_row","array_col"]].max().to_dict())
print(open(f"{B}/square_002um/spatial/scalefactors_json.json").read())
