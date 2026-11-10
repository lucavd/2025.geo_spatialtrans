#!/usr/bin/env python3
# tools/S2.7_mutants.py - S2.7, G-1: otto mutanti della simulazione a pezzi; ognuno deve far fallire la batteria R/testing/test_S2.7.R
# (gruppo a su linfonodo, cervello, cuore e dischi; il cervello ha il fondo). Uso (root del repo): .venv/bin/python tools/S2.7_mutants.py [processi]
import os, shutil, subprocess, sys, tempfile, concurrent.futures as cf
F6 = "05b6_simulate_visium_hd_pieces.R"
MUT = [
 ("M1 flussi per pezzo invece che per tessera", F6, "assign(\".Random.seed\", env$st_place[[ut[h]]], envir = globalenv())", "assign(\".Random.seed\", env$st_place[[k]], envir = globalenv())"),
 ("M2 trascritti in uscita scartati", F6, "for (d in sort(unique(c(dt, do_)))) saveRDS(", "for (d in intersect(k, sort(unique(c(dt, do_))))) saveRDS("),
 ("M3 alone H = 5 D", F6, "env$halo_um <- 25 * p$D", "env$halo_um <- 5 * p$D"),
 ("M4 quadratino di confine in due pezzi", F6, "c(i0 = pr * s, i1 = min(grid$n_row, (pr + 1L) * s) - 1L,", "c(i0 = pr * s, i1 = min(grid$n_row - 1L, (pr + 1L) * s),"),
 ("M5 campo senza alone (rumore nullo fuori dalla finestra)", F6, "kr <- which(gr >= r0 & gr < r0 + NR); kc <- which(gc >= c0 & gc < c0 + NC)", "kr <- which(gr >= P0 & gr < P0 + nr); kc <- which(gc >= Q0 & gc < Q0 + nc)"),
 ("M6 tetto dal massimo del pezzo", F6, "kt <- stats::rbinom(length(it), tr$n[it], ct[tr$bin_id[it]] / Ms)", "kt <- stats::rbinom(length(it), tr$n[it], ct[tr$bin_id[it]] / max(ct[ids]))"),
 ("M7 composizione del fondo per pezzo", F6, "gb <- sample.int(G, sum(K), replace = TRUE, prob = env$comp)", "gb <- sample.int(G, sum(K), replace = TRUE, prob = tabulate(tr$gene, G) + 1)"),
 ("M8 seconda coda ripiegata sul riquadro della tessera", F6, "if (env$frame_open) { fb <- env$frame; x1[sec]", "if (env$frame_open) { fb <- c(xmin = min(x), ymin = min(y), xmax = max(x) + 1e-9, ymax = max(y) + 1e-9); x1[sec]"),
]
def one(m):
    name, f, old, new = m
    d = tempfile.mkdtemp(prefix="s27mut_")
    for x in os.listdir("R"):
        if x.endswith(".R"): shutil.copy(os.path.join("R", x), d)
    p = os.path.join(d, f); s = open(p, encoding="utf-8").read()
    assert s.count(old) == 1, (name, s.count(old))
    open(p, "w", encoding="utf-8").write(s.replace(old, new))
    env = dict(os.environ, S27_RDIR=d, S27_TISSUES="A4,A5,A6", R_LIBS_SITE="/nonexistent", R_LIBS_USER="/nonexistent", OMP_NUM_THREADS="1")
    r = subprocess.run(["Rscript", "--vanilla", "R/testing/test_S2.7.R", "a"], env=env, capture_output=True, text=True)
    fails = [l for l in r.stdout.splitlines() if l.startswith("FAIL")]
    shutil.rmtree(d, ignore_errors=True)
    return name, r.returncode != 0, fails[:3], r.stderr[-300:] if r.returncode != 0 and not fails else ""
w = int(sys.argv[1]) if len(sys.argv) > 1 else 4
with cf.ThreadPoolExecutor(w) as ex: res = list(ex.map(one, MUT))
nd = 0
for name, det, fails, err in res:
    nd += det
    print(("RILEVATO " if det else "NON RILEVATO ") + name + (" | " + " ; ".join(fails) if fails else "") + (" | errore: " + err.replace("\n", " ") if err else ""))
print(f"mutanti rilevati: {nd}/{len(MUT)}")
sys.exit(0 if nd == len(MUT) else 1)
