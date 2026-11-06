"""tools/S2.4_mutants.py — S2.4, G-1: mutanti del codice (copie modificate di R/ passate alla batteria con S24_RSRC).
Ogni mutante deve far fallire almeno un test. Uso: .venv/bin/python tools/S2.4_mutants.py <nproc> -> results/S2.4/S2.4_mutants.txt"""
import os, sys, shutil, subprocess, tempfile
from multiprocessing import Pool
F4, F3 = "R/05b4_capture_and_bin.R", "R/05b3_place_and_move_transcripts.R"
M = [
 ("assottigliamento con c invece di c/C", F4, "p <- cvec / C", "p <- pmin(1, cvec)"),
 ("normalizzazione su tutti i quadratini", F4, "nrm <- mean(craw[tissue])", "nrm <- mean(craw)"),
 ("parita' sulle righe pari", F4, "odd_r <- (seq_len(nr) - 1L) %% 2L == 1L", "odd_r <- (seq_len(nr) - 1L) %% 2L == 0L"),
 ("segno di theta", F4, "th <- if (rho == 0) 0 else (1 - sqrt(max(0, 1 - 4 * rho^2))) / (2 * rho)", "th <- if (rho == 0) 0 else -(1 - sqrt(max(0, 1 - 4 * rho^2))) / (2 * rho)"),
 ("fondo senza strisce", F4, "K <- stats::rpois(nb, params$b * cvec)", "K <- stats::rpois(nb, params$b + 0 * cvec)"),
 ("raggio della seconda coda esponenziale", F3, "rho2 <- stats::rgamma(nS, shape = 2, scale = params$ell2)", "rho2 <- stats::rgamma(nS, shape = 1, scale = params$ell2)"),
 ("riflessione mancante", F3, "if (frame_open) { x1[sec] <- .s24_fold", "if (FALSE) { x1[sec] <- .s24_fold"),
 ("seconda coda solo fuori dalla prima (miscela sbagliata)", F3, "sec <- stats::runif(nt) < w2; nS <- sum(sec)", "sec <- stats::runif(nt) < w2 & mode != 2L; nS <- sum(sec)"),
 ("blocchi di 8 um sfasati di un quadratino", F4, "X8 <- sum_to_8um(X2, g)", "X8 <- sum_to_8um(X2[, c(2:nb, 1)], g)"),
 ("seme ignorato nei blocchi", F4, 'assign(".Random.seed", ss[[k]], envir = globalenv())', 'assign(".Random.seed", ss[[1]], envir = globalenv())'),
 ("ipergeometrica con nucleo e citoplasma scambiati", F4, "kn <- stats::rhyper(length(io), og$n_born_nucleus[io], og$n[io] - og$n_born_nucleus[io], kk)", "kn <- stats::rhyper(length(io), og$n[io] - og$n_born_nucleus[io], og$n_born_nucleus[io], kk)"),
 ("seme non rimescolato", F4, 'RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); s1 <- sample.int(.Machine$integer.max, 1L)', "s1 <- seed"),
]
def run(k):
    name, f, a, b = M[k]; d = tempfile.mkdtemp(prefix="s24mut_"); shutil.copytree("R", f"{d}/R")
    s = open(f"{d}/{f}").read(); assert s.count(a) == 1, (name, s.count(a)); open(f"{d}/{f}", "w").write(s.replace(a, b))
    env = dict(os.environ, S24_RSRC=d, R_LIBS_SITE="/nonexistent", R_LIBS_USER="/nonexistent", OMP_NUM_THREADS="1")
    p = subprocess.run(["Rscript", "--vanilla", "R/testing/test_S2.4.R"], capture_output=True, text=True, env=env)
    fails = [l for l in p.stdout.splitlines() if l.startswith("FAIL")]
    shutil.rmtree(d)
    return name, len(fails), [l.split(" ")[1] for l in fails]
if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    with Pool(n) as p: R = p.map(run, range(len(M)), chunksize=1)
    out = [f"{'RILEVATO' if k > 0 else 'NON RILEVATO'} | {nm} | test falliti: {k} ({', '.join(f)})" for nm, k, f in R]
    det = sum(k > 0 for _, k, _ in R); out.append(f"mutanti rilevati {det}/{len(R)}")
    open("results/S2.4/S2.4_mutants.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
