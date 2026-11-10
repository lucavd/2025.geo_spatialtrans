"""tools/S2.6_mutants.py - S2.6, G-1: mutanti degli strumenti della prova di prestazione; ognuno deve far fallire R/testing/test_S2.6.R.
Uso (root del repo): .venv/bin/python tools/S2.6_mutants.py [out.txt]"""
import os, sys, shutil, subprocess, tempfile
P, M = "S2.6_pieces.R", "S2.6_measures.R"
MUT = [
    ("ancore tolte (griglia del pezzo non allineata dove il tessuto non arriva al riquadro)", P, "P$region_polygons <- c(Rp, anc)", "P$region_polygons <- Rp"),
    ("taglio dei territori al riquadro saltato", P, "  if (length(cross)) {\n    Tt[cross]", "  if (FALSE) {\n    Tt[cross]"),
    ("quadratini tenuti con una colonna in meno", P, "kj <- jj[xc >= U[1] & xc < U[3]]", "kj <- jj[xc > U[1] + 2 & xc < U[3]]"),
    ("ogni trascritto attribuito al pezzo proprietario", P, "og$own <- og$cell_id %in% own_ids", "og$own <- TRUE"),
    ("stesso seme in tutti i pezzi", P, "sd_p <- s26_seeds(seed, nrow(P))", "sd_p <- rep(s26_seeds(seed, 1L), nrow(P))"),
    ("blocchi di 8 um sfasati nella cucitura", P, "Sm <- Matrix::sparseMatrix(i = seq_along(gid), j = b8,", "Sm <- Matrix::sparseMatrix(i = seq_along(gid), j = pmin(n8, b8 + 1L),"),
    ("riflessione spenta nei pezzi", P, "frame_open = TRUE, layers = s26_layers(st), workers = 1L)", "frame_open = FALSE, layers = s26_layers(st), workers = 1L)"),
    ("giunture spostate di 8 um nelle misure", M, "sx <- xmin + ks * Tu; sy <- ymin + ks * Tu", "sx <- xmin + ks * Tu + 8; sy <- ymin + ks * Tu + 8"),
]
out = sys.argv[1] if len(sys.argv) > 1 else "results/S2.6/S2.6_mutants.txt"
os.makedirs(os.path.dirname(out), exist_ok=True)
env = dict(os.environ, R_LIBS_SITE="/nonexistent", R_LIBS_USER="/nonexistent", OMP_NUM_THREADS="1")
lines = []; ok_all = True
for name, f, a, b in MUT:
    src = open(f"tools/{f}").read(); n = src.count(a)
    if n != 1: lines.append(f"ERRORE {name}: frammento trovato {n} volte"); ok_all = False; continue
    d = tempfile.mkdtemp(prefix="s26mut_")
    for g in (P, M, "S2.6_section_layer.R"): shutil.copy(f"tools/{g}", f"{d}/{g}")
    open(f"{d}/{f}", "w").write(src.replace(a, b))
    p = subprocess.run(["Rscript", "--vanilla", "R/testing/test_S2.6.R"], env=dict(env, S26_TOOLS=d), capture_output=True, text=True)
    fails = [l for l in p.stdout.splitlines() if l.startswith("FAIL")]
    det = p.returncode != 0
    if not fails and det: fails = ["errore R: " + (p.stderr.strip().splitlines() or ["?"])[-1]]
    ok_all &= det
    lines.append(f"{'RILEVATO' if det else 'NON RILEVATO'} {name}: {len(fails)} asserzioni fallite; {'; '.join(x[:90] for x in fails[:3])}")
    shutil.rmtree(d)
lines.append(f"TOTALE: {sum(l.startswith('RILEVATO') for l in lines)}/{len(MUT)} rilevati")
open(out, "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
sys.exit(0 if ok_all else 1)
