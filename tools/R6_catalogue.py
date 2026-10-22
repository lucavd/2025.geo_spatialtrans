"""R6 catalogue (same content as in tools/R6_cluster.py; shared by R6_relabel.py)."""
SPECIES = {"A1": "mouse", "A2": "human", "A3": "human", "A4": "human", "A5": "mouse", "A6": "human"}
# coarse catalogue (approved by Luca 2026-10-23) -> PanglaoDB cell types
CAT = {
 "A1": {"epiteliale": ["Enterocytes", "Goblet cells", "Paneth cells", "Enteroendocrine cells", "Tuft cells", "Crypt cells", "Epithelial cells"],
        "stromale e muscolo liscio": ["Fibroblasts", "Smooth muscle cells", "Stromal cells", "Myofibroblasts", "Pericytes"],
        "immunitario": ["B cells", "T cells", "Plasma cells", "Macrophages", "Dendritic cells", "Mast cells", "Monocytes", "NK cells", "Neutrophils", "T memory cells", "T helper cells"],
        "endoteliale": ["Endothelial cells"]},
 "A2": {"epiteliale": ["Enterocytes", "Goblet cells", "Crypt cells", "Epithelial cells", "Enteroendocrine cells", "Tuft cells"],
        "fibroblasti": ["Fibroblasts", "Myofibroblasts", "Stromal cells"],
        "immunitario": ["B cells", "T cells", "Plasma cells", "Macrophages", "Dendritic cells", "Mast cells", "Monocytes", "NK cells", "Neutrophils", "T memory cells", "T helper cells"],
        "endoteliale": ["Endothelial cells"],
        "muscolo liscio": ["Smooth muscle cells", "Pericytes"]},
 "A4": {"linfociti B": ["B cells", "B cells naive", "B cells memory"],
        "linfociti T": ["T cells", "T helper cells", "T memory cells", "T cytotoxic cells", "T regulatory cells", "T follicular helper cells"],
        "plasmacellule": ["Plasma cells"],
        "mieloidi": ["Macrophages", "Dendritic cells", "Monocytes", "Plasmacytoid dendritic cells", "Mast cells", "Neutrophils"],
        "stromali": ["Fibroblasts", "Stromal cells", "Smooth muscle cells", "Pericytes"],
        "endoteliali": ["Endothelial cells"]},
 "A5": {"neuroni": ["Neurons", "Interneurons", "Pyramidal cells", "GABAergic neurons", "Glutaminergic neurons"],
        "astrociti": ["Astrocytes"], "oligodendrociti": ["Oligodendrocytes", "Oligodendrocyte progenitor cells"], "microglia": ["Microglia"],
        "vascolari": ["Endothelial cells", "Endothelial cells (blood brain barrier)", "Pericytes", "Smooth muscle cells"]},
 "A6": {"cardiomiociti": ["Cardiomyocytes"], "fibroblasti": ["Fibroblasts"], "endoteliali": ["Endothelial cells", "Endothelial cells (aorta)"],
        "immunitari": ["Macrophages", "T cells", "B cells", "Mast cells", "Monocytes", "Dendritic cells", "Neutrophils"],
        "muscolo liscio e periciti": ["Smooth muscle cells", "Vascular smooth muscle cells", "Pericytes"]}}
CAT["A3"] = CAT["A2"]
# tumour vs normal epithelium (A2/A3 only): crypt/stem-like score vs differentiated (enterocyte + goblet) score; technical rule, declared
EPI_SPLIT = {"tumorale": ["Crypt cells"], "normale": ["Enterocytes", "Goblet cells"]}
