
"""R5 exploration E0 (geometry only, no counts): how much tissue edge facing glass each section has.
Reads only tissue_positions.parquet (in_tissue) of the 2 um grid. Declared before pre-registration."""
import sys, json, numpy as np, pandas as pd
from scipy import ndimage as ndi
DS = ["Visium_HD_Mouse_Small_Intestine","Visium_HD_Human_Colon_Cancer","Visium_HD_Human_Lymph_Node_FFPE",
      "Visium_HD_6p5mm_Mouse_Brain","Visium_HD_6p5mm_Human_Heart"]
out = []
for ds in DS:
    p = f"data/real/datasets/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet"
    t = pd.read_parquet(p, columns=["array_row","array_col","in_tissue"])
    nr, nc = t.array_row.max()+1, t.array_col.max()+1
    m = np.zeros((nr,nc), bool); m[t.array_row.values, t.array_col.values] = t.in_tissue.values.astype(bool)
    # distance (um) of out-of-tissue bins to tissue
    d_out = ndi.distance_transform_edt(~m) * 2.0
    lab, nlab = ndi.label(~m)
    border_labels = set(np.unique(np.r_[lab[0,:],lab[-1,:],lab[:,0],lab[:,-1]])) - {0}
    holes = [l for l in range(1, nlab+1) if l not in border_labels]
    sizes = ndi.sum(np.ones_like(lab), lab, index=holes) if holes else np.array([])
    # boundary bins of tissue facing out-of-tissue (4-neighbour)
    er = ndi.binary_erosion(m, structure=np.array([[0,1,0],[1,1,1],[0,1,0]]), border_value=1)
    edge = m & ~er
    # edge bins whose outward glass band is wide: max out-distance within 160 um square around exceeds 150 um
    far = ndi.maximum_filter(d_out, size=121)  # +-60 bins = +-120 um window
    edge_wide = edge & (far >= 100)
    # array border touching tissue
    border = np.r_[m[0,:], m[-1,:], m[:,0], m[:,-1]]
    out.append(dict(dataset=ds, grid=f"{nr}x{nc}", frac_in_tissue=round(m.mean(),4),
        tissue_mm2=round(m.sum()*4e-6,2), glass_ge50um_mm2=round((d_out>=50).sum()*4e-6,2),
        edge_len_mm=round(edge.sum()*2e-3,1), edge_len_glass100_mm=round(edge_wide.sum()*2e-3,1),
        frac_array_border_in_tissue=round(border.mean(),3),
        n_holes_sr=len(holes), n_holes_ge100um2=int((sizes*4>=100).sum()), n_holes_ge2500um2=int((sizes*4>=2500).sum())))
df = pd.DataFrame(out); df.to_csv(sys.argv[1], index=False); print(df.to_string())
