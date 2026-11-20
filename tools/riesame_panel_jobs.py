"""tools/riesame_panel_jobs.py - riesame, Parte B: liste delle esecuzioni della catena (tools/riesame_panel_run.sh) per i bracci V, T, S, B
e per il cancello G-4. Aree: V = 5 aree di taratura P2 per tumore (A2) e stroma (A3) (F3.1/main/real); T = aree L2 di F3.3 del tumore
e dello stroma (P1, P5, tumore da 6,5 mm; righe OUT/L2 di results/F3.3/F3.3_area_values.csv, esportazioni F3.3/runs/<run>);
S = simulazioni di F3.1 (F3.1/main/sim, 5 aree x 2 semi); B = simulazioni banali (riesame/B/sim, tools/riesame_B_sim.R).
G-4: un'area per sezione rifatta con la catena sui geni del confronto di F3.3 (V e T: common; S: all, i 18 132 geni con cui F3.3 ha misurato S).
Uso: python3 tools/riesame_panel_jobs.py  ->  /mnt/micron/geo_spatialtrans/riesame/B/jobs/*.tsv"""
import os, csv
M = "/mnt/micron/geo_spatialtrans"; OUT = f"{M}/riesame/B/jobs"; os.makedirs(OUT, exist_ok=True)
AV = list(csv.DictReader(open(os.path.expanduser("~/2025.geo_spatialtrans/results/F3.3/F3.3_area_values.csv"))))
def V(genes, arm="V", areas=None):
    return [(arm, genes, f"{M}/F3.1/main/real/{A}_r{k}", A, "real", f"r{k}", f"{A}_r{k}", f"{A}_r{k}", "IN", f"{A}_IN", "IN", f"{A}_r{k}")
            for A in ("A2", "A3") for k in range(1, 6) if areas is None or f"{A}_r{k}" in areas]
def T(genes, arm="T", runs=None):
    return [(arm, genes, f"{M}/F3.3/runs/{r['run']}", r["A"], "real", "-", r["unit"], r["run"], "OUT", r["set"], r["level"], r["unit"])
            for r in AV if r["src"] == "OUT" and r["level"] == "L2" and r["A"] in ("A2", "A3") and (runs is None or r["run"] in runs)]
def S(genes, arm="S", runs=None, root=f"{M}/F3.1/main/sim", src="SIM"):
    return [(arm, genes, f"{root}/{A}_r{k}_k{s}", A, "sim", "-", f"{A}_r{k}", f"{A}_r{k}_k{s}", src, f"{A}_{src}", src, f"{A}_r{k}")
            for A in ("A2", "A3") for k in range(1, 6) for s in (0, 1) if runs is None or f"{A}_r{k}_k{s}" in runs]
def write(name, rows):
    with open(f"{OUT}/{name}.tsv", "w") as f:
        for r in rows: f.write("\t".join(r) + "\n")
    print(name, len(rows))
G4T = ["CRC_P1_SR300_tumore_2_cg", "CRC_P5_SR300_tumore_1_cg", "CRC_6p5mm_SR410_tumore_2_cg"]
write("g4", V("common", "G4", {"A2_r1", "A3_r1"}) + T("common", "G4", set(G4T)) + S("all", "G4", {"A2_r1_k0", "A3_r1_k0"}))
write("panel_V", V("panel")); write("panel_T", T("panel")); write("panel_S", S("panel")); write("common_S", S("common"))
write("panel_B", S("panel", "B", root=f"{M}/riesame/B/sim", src="B")); write("common_B", S("common", "B", root=f"{M}/riesame/B/sim", src="B"))
