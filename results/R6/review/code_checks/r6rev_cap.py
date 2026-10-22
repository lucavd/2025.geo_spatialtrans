
import pandas as pd, os, sys
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0,"tools"); from R6_catalogue import CAT
M="/mnt/micron/geo_spatialtrans/R6"
pan = pd.read_csv(f"{M}/PanglaoDB_omnipath_mouse.tsv", sep="\t"); pan = pan[pan.mouse==True]
sym = pan.genesymbol.astype(str); ch = sym != sym.str.capitalize()
print("mouse rows", len(pan), "symbols changed by capitalize():", int(ch.sum()), "unique:", sym[ch].nunique(), "examples:", sorted(set(sym[ch]))[:25])
for A, key in [("A1","A1"),("A5","A5")]:
    genes = set(pd.read_csv(f"{M}/counts/{key}_genes.tsv", sep="\t").gene_name)
    for ty, pts in CAT[A].items():
        s = pan.genesymbol[pan.cell_type.isin(pts)].astype(str)
        lost = sorted(set(s[(s != s.str.capitalize()) & s.isin(genes) & ~s.str.capitalize().isin(genes)]))
        print(f"  {A} {ty}: marker symbols present in h5 as written but lost after capitalize(): {len(lost)} {lost[:12]}")
