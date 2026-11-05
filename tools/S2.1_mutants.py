"""tools/S2.1_mutants.py - S2.1, cancello G-1 (e): 8 mutanti di R/05b1_build_hd_grid.R, ognuno deve far fallire la batteria.
Uso (root del repo): python3 -B tools/S2.1_mutants.py <outcsv>. Ogni mutante e' una copia di R/ in una cartella temporanea con UNA
sostituzione di testo (verificata: deve comparire esattamente una volta); la batteria gira con S21_RSRC=<cartella>."""
import csv, os, shutil, subprocess, sys, tempfile
SRC = "R/05b1_build_hd_grid.R"
MUT = [
    ("M1", "i e j scambiati nei centri", "x = grid$x0 + (j + 0.5) * grid$bin_um, y = grid$y0 + (i + 0.5) * grid$bin_um",
     "x = grid$x0 + (i + 0.5) * grid$bin_um, y = grid$y0 + (j + 0.5) * grid$bin_um"),
    ("M2", "centro spostato di mezzo quadratino", "x = grid$x0 + (j + 0.5) * grid$bin_um", "x = grid$x0 + (j + 1) * grid$bin_um"),
    ("M3", "margine non aggiunto", "x0 <- bb[1] - margin_um; y0 <- bb[2] - margin_um", "x0 <- bb[1]; y0 <- bb[2]"),
    ("M4", "quadratini sul bordo di una tessera contati due volte", "i1 = pmin(ri[g$h] + nt, grid$n_row)", "i1 = pmin(ri[g$h] + nt + 1L, grid$n_row)"),
    ("M5", "blocchi di 3 invece di 4 nella tabella dei quadratini (revisione F05)", "i8 <- i %/% grid$block; j8 <- j %/% grid$block", "i8 <- i %/% 3L; j8 <- j %/% 3L"),
    ("M6", "nucleo calcolato sul territorio", "nuclei[cc[keep]], grid_um = grid_um", "territories[cc[keep]], grid_um = grid_um"),
    ("M7", "in_tissue con > invece di >=", "bins$tissue_frac >= 0.5", "bins$tissue_frac > 0.5"),
    ("M8", "ultima riga e colonna perse (arrotondamento per difetto)", "max(block, block * ceiling(n / block))", "max(block, block * floor(n / block))"),
    ("M9", "coppie contate in due tessere della verita' geometrica (revisione F06)", "k <- which(tile == t); i <- ci[k]", "k <- which(tile == t | tile == t - 1L); i <- ci[k]"),
    ("M10", "cell_id della verita' non riportato a quello dello Step 1 (revisione F07)", "  ov$cell_id <- cell_layer$cell_df$cell_id[ov$cell_id]\n", "  ov$cell_id <- ov$cell_id\n"),
    ("M11", "soglia di in_tissue a 0.45 (revisione F08)", "bins$tissue_frac >= 0.5", "bins$tissue_frac >= 0.45"),
]
# Non rilevabile per costruzione, dichiarato (revisione F04): l'origine globale invece di quella della tessera cambia le aree di ~1e-12 um2 a 1e4 um,
# sotto ogni tolleranza; BL-104 si chiude con la misura del revisore, non con un mutante.
out = sys.argv[1] if len(sys.argv) > 1 else "results/S2.1/S2.1_mutants.csv"
os.makedirs(os.path.dirname(out), exist_ok=True)
code = open(SRC).read(); rows = []
env = dict(os.environ, R_LIBS_SITE="/nonexistent", R_LIBS_USER="/nonexistent")
for mid, what, a, b in MUT:
    n = code.count(a)
    if n != 1:
        rows.append(dict(mutant=mid, what=what, applied=False, n_fail=None, crashed=None, detected=False)); continue
    d = tempfile.mkdtemp(prefix=f"s21_{mid}_"); os.makedirs(f"{d}/R")
    for f in os.listdir("R"):
        if f.endswith(".R"): shutil.copy(f"R/{f}", f"{d}/R/{f}")
    open(f"{d}/R/05b1_build_hd_grid.R", "w").write(code.replace(a, b))
    p = subprocess.run(["Rscript", "--vanilla", "R/testing/test_S2.1.R"], env=dict(env, S21_RSRC=d), capture_output=True, text=True)
    fails = [l for l in p.stdout.splitlines() if l.startswith("FAIL")]
    crashed = p.returncode != 0 and not fails
    rows.append(dict(mutant=mid, what=what, applied=True, n_fail=len(fails), crashed=crashed, detected=bool(fails),
                     first_fail=(fails[0][:90] if fails else (p.stderr.strip().splitlines() or [""])[-1][:90])))
    shutil.rmtree(d)
    print(mid, rows[-1]["detected"], rows[-1]["n_fail"], rows[-1]["first_fail"], flush=True)
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["mutant", "what", "applied", "n_fail", "crashed", "detected", "first_fail"]); w.writeheader(); w.writerows(rows)
print(f"mutanti rilevati da un FAIL (i crash non contano): {sum(r['detected'] for r in rows)}/{len(MUT)}; crash: {sum(bool(r['crashed']) for r in rows)}")
