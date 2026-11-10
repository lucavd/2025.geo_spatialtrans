# 05_deconvolve.R — CPU deconvolution of the main datasets written by 04_simulate.R (PROTOCOLLO §6 T5;
# ADEMP §4, §5.2, §7.2, §8). cell2location (GPU, server) reads the same files with its own script.
#
# Design (fixed by ADEMP v1.2; implementation written 2026-11-05 before the first run):
#   one job = one replicate x one condition (6 datasets: anchor, E, C−, E′, O, P) x 4 methods
#   (RCTD full, spatialDWLS, MuSiC, NNLS; ADEMP §4 settings, R/deconvolve.R).
#   Reference: the real capped training reference for every dataset except P; for P the reference is
#   rebuilt from world_p.rds (reference_p_from_world(), R/arm_p.R), with NNLS profiles = mean counts
#   per simulated reference cell (ADEMP §2.1).
#   Outputs cleaned with clean_props() (ADEMP §7.2: tolerance 1e-5, values in [-1e-5, 0) set to 0,
#   rows rescaled); failed spots get the null-baseline proportions (share of each type among all
#   training cells) in the main analysis; truth = cell shares (main), UMI shares (sensitivity (c)).
#   Seeds: each method run starts from its own sub-stream of the replicate stream of 04_simulate.R,
#   numbered 1000 + 100 × condition + 10 × dataset slot + method (disjoint from the generation
#   sub-streams, which stop at 272), so a dataset or method added later changes nothing else.
#
# Pre-registered checks (written 2026-11-05 before the first run; the job stops if one fails):
#   D1. spots of every dataset have the genes of their reference (P: P00001…; others: real genes);
#   D2. the rebuilt P reference has the cells, types and donors of the real reference;
#   D3. every method returns (after completion) 500 spots x all reference types; after cleaning, the
#       returned rows have no negative value and sum to 1 (1e-12);
#   D4. RMSE and JSD are computed on the filled matrix (no NA).
# Reported per job: RMSE (filled, cell truth), JSD, failure rate, RMSE on the UMI truth, RMSE on the
# spots returned by all methods (sensitivity (a)), seconds, error message, smallest raw value
# (residue check of ADEMP §7.2), RNG state after the method (seed-reset check, ADEMP §8).
#
# Usage (working directory = thesis/):
#   Rscript scripts/05_deconvolve.R <tissue> [pilot|main] [n_rep] [workers]
# Output: <TESI_DATA>/sim/<mode>/<tissue>/r<rep>/est_c<cells>_d<depth>.rds and <tissue>_deconv_metrics.csv.

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
ref_mode <- if (length(args) >= 2) args[2] else "pilot"
n_rep <- if (length(args) >= 3) as.integer(args[3]) else 81L
workers <- if (length(args) >= 4) as.integer(args[4]) else 8L
stopifnot(tissue %in% c("colon", "lymphoid", "brain"), ref_mode %in% c("pilot", "main"))

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
src <- c("R/simulate_spots.R", "R/arm_p.R", "R/deconvolve.R", "R/pilot_metrics.R")
suppressPackageStartupMessages(library(Matrix)); for (f in src) source(f)
R0 <- readRDS(file.path(data_dir, "derived", "reference", sprintf("%s_%s.rds", tissue, ref_mode)))
sim_dir <- file.path(data_dir, "sim", ref_mode, tissue)
log_file <- file.path(dirname(sim_dir), paste0(tissue, "_deconv_progress.log"))
methods <- c("NNLS", "RCTD", "MuSiC", "spatialDWLS")
slot <- c(anchor = 1, E = 2, Cminus = 3, Eprime = 4, O = 5, P = 6, G = 7)
types <- R0$types_ref
p_null <- null_props(lengths(R0$L_by_type))[types]
p_null_umi <- vapply(R0$L_by_type, sum, numeric(1)); p_null_umi <- (p_null_umi / sum(p_null_umi))[types]
L_ref <- Matrix::colSums(R0$ref)
cond_keys <- c("c02_d0175", "c02_d4900", "c10_d0175", "c10_d4900")
jobs <- expand.grid(cond = cond_keys, rep = seq_len(n_rep), stringsAsFactors = FALSE)

substream <- function(seed, k) { for (j in seq_len(k)) seed <- parallel::nextRNGSubStream(seed); seed }

run_job <- function(r) {
  key <- jobs$cond[r]; i <- jobs$rep[r]
  rd <- file.path(sim_dir, sprintf("r%03d", i))
  f_out <- file.path(rd, paste0("est_", key, ".rds"))
  if (file.exists(f_out)) return(readRDS(f_out)$metrics)
  log_line <- function(msg) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = log_file, append = TRUE)
  log_line(sprintf("start rep %d %s", i, key))
  x <- readRDS(file.path(rd, paste0(key, ".rds")))
  sp <- rownames(x$truth); nc <- x$meta$n_cells
  refs <- list(real = list(ref = R0$ref, meta = R0$ref_meta, profiles = R0$profiles))
  if ("P" %in% names(x$counts)) {
    RP <- reference_p_from_world(readRDS(file.path(rd, "world_p.rds")), R0$ref_meta, L_ref, colnames(R0$ref))
    if (ncol(RP$ref) != ncol(R0$ref) || !identical(RP$ref_meta, R0$ref_meta)) stop("D2")
    refs$P <- list(ref = RP$ref, meta = RP$ref_meta, profiles = RP$profiles)
  }
  ci <- match(key, cond_keys)
  rows <- list(); est_all <- list()
  for (d in names(x$counts)) {
    rf <- refs[[if (d == "P") "P" else "real"]]
    spots <- x$counts[[d]]
    if (!identical(rownames(spots), rownames(rf$ref))) stop("D1: ", d)
    est <- list()
    for (mi in seq_along(methods)) {
      m <- methods[mi]
      assign(".Random.seed", substream(x$meta$rep_seed, 1000 + 100 * ci + 10 * slot[[d]] + mi), envir = .GlobalEnv)
      out <- tryCatch(switch(m,
        NNLS = run_nnls(spots, rf$profiles),
        RCTD = run_rctd(spots, rf$ref, rf$meta$cell_type),
        MuSiC = run_music(spots, rf$ref, rf$meta$cell_type, rf$meta$donor_id),
        spatialDWLS = run_dwls(spots, rf$ref, rf$meta$cell_type, n_cell = nc)),
        error = function(e) list(prop = NULL, seconds = NA, error = conditionMessage(e)))
      p <- if (is.null(out$prop)) matrix(NA_real_, length(sp), length(types), dimnames = list(sp, types))
           else complete_prop(out$prop, sp, types)
      p <- clean_props(p)
      if (!identical(dim(p), c(length(sp), length(types)))) stop("D3 dims: ", d, " ", m)
      okr <- !failed_spots(p)
      if (any(okr) && (any(p[okr, ] < 0) || any(abs(rowSums(p[okr, , drop = FALSE]) - 1) > 1e-12))) stop("D3 rows: ", d, " ", m)
      est[[m]] <- p
      rows[[length(rows) + 1]] <- data.frame(dataset = d, method = m, seconds = out$seconds,
        error = if (is.null(out$error)) NA else out$error, fail_rate = mean(!okr),
        min_raw = attr(p, "min_raw"), rng_after = sum(as.numeric(.Random.seed)), stringsAsFactors = FALSE)
    }
    fs <- lapply(est, failed_spots)
    fs <- fs[!vapply(fs, all, logical(1))]
    common <- if (length(fs)) !Reduce(`|`, fs) else rep(FALSE, length(sp))   # sensitivity (a)
    tu <- x$truth_umi[[d]][, types]
    for (m in methods) {
      k <- which(vapply(rows, function(z) z$dataset == d && z$method == m, logical(1)))
      fl <- fill_failed(est[[m]], p_null)
      if (anyNA(fl)) stop("D4: ", d, " ", m)
      rows[[k]]$rmse <- rmse(fl, x$truth); rows[[k]]$jsd <- jsd(fl, x$truth)
      rows[[k]]$rmse_umi <- rmse(fill_failed(est[[m]], p_null_umi), tu)
      rows[[k]]$n_common <- sum(common)
      rows[[k]]$rmse_common <- if (any(common)) rmse(est[[m]][common, ], x$truth[common, ]) else NA
    }
    nb <- matrix(p_null, length(sp), length(types), byrow = TRUE, dimnames = list(sp, types))
    rows[[length(rows) + 1]] <- data.frame(dataset = d, method = "null baseline", seconds = 0, error = NA,
      fail_rate = 0, min_raw = NA, rng_after = NA, rmse = rmse(nb, x$truth), jsd = jsd(nb, x$truth),
      rmse_umi = rmse(matrix(p_null_umi, length(sp), length(types), byrow = TRUE), tu), n_common = NA, rmse_common = NA)
    est_all[[d]] <- est
  }
  met <- cbind(tissue = tissue, rep = i, n_cells = nc, depth = x$meta$depth, do.call(rbind, rows))
  tmp <- paste0(f_out, ".tmp"); saveRDS(list(metrics = met, props = est_all), tmp); file.rename(tmp, f_out)
  log_line(sprintf("done rep %d %s", i, key))
  met
}

cat(sprintf("[%s, %s] %d jobs (%d replicates x 4 conditions), %d workers\n", tissue, ref_mode, nrow(jobs), n_rep, workers))
t0 <- proc.time()[["elapsed"]]
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")   # one core per worker (ADEMP §6); workers inherit it
cl <- parallel::makeCluster(workers, setup_strategy = "sequential")
lib <- .libPaths()
invisible(parallel::clusterCall(cl, function(lib) { .libPaths(lib); RNGkind("L'Ecuyer-CMRG"); NULL }, lib))
parallel::clusterExport(cl, c("R0", "sim_dir", "log_file", "methods", "slot", "types", "p_null", "p_null_umi", "L_ref",
                              "cond_keys", "jobs", "substream", "run_job", "tissue", "src"))
invisible(parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(Matrix)); for (f in src) source(f); NULL }))
done <- parallel::parLapplyLB(cl, seq_len(nrow(jobs)), function(r)
  tryCatch(run_job(r), error = function(e) paste("ERROR job", r, conditionMessage(e))))
parallel::stopCluster(cl)
errs <- unlist(Filter(is.character, done))
if (length(errs)) { cat(errs, sep = "\n"); stop(length(errs), " jobs failed a check") }
m <- do.call(rbind, done)
write.csv(m, file.path(dirname(sim_dir), paste0(tissue, "_deconv_metrics.csv")), row.names = FALSE)
cat(sprintf("Done in %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
cat("Mean RMSE (cell truth, failed spots filled) by dataset x method:\n")
print(round(tapply(m$rmse, list(m$dataset, m$method), mean), 3))
cat("Failure rate:\n"); print(round(tapply(m$fail_rate, list(m$dataset, m$method), mean), 3))
cat(sprintf("Smallest raw value over all runs: %.2e (tolerance 1e-5)\n", min(m$min_raw, na.rm = TRUE)))
if (any(duplicated(m$rng_after[!is.na(m$rng_after)]))) cat("WARNING: identical RNG state after two method runs (seed reset?)\n")
