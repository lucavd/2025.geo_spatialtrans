# 06_pilot.R — declared pilot (PROTOCOLLO §6 week 2; ADEMP §10). Its numbers never enter the
# thesis results: it only fills the [PILOT] placeholders (run time, failures, SD(RMSE_i) for n_sim,
# seed check, reuse of held-out cells).
#
# Design (ADEMP §10, proposal 2026-10-29; levels Pierangela 2026-10-30):
#   pilot donors from 03_reference.R <tissue> pilot (drawn inside training, seed 20261030);
#   arms E and C− (one new gene permutation per replicate), plus the real anchor (built from the
#   pilot held-out donor), because n_sim takes the largest SD across the design (ADEMP §6) and
#   the anchor is part of every condition;
#   cells per spot 2 / 10 (synthspot bounds; 2 = Visium HD 8 um bin, 10 = standard Visium);
#   depth per cell 175 / 4,900 UMI (10x Visium HD and standard Visium colorectal datasets, Space
#   Ranger metrics); spot target = cells x UMI per cell, spots above are thinned (Pierangela 2026-10-31);
#   composition Dirichlet(alpha_k = K * p_k), p_k = share of spot type k among training cells, then
#   cells drawn multinomially (Pierangela 2026-10-31);
#   500 spots per dataset; 10 replicates per condition;
#   methods RCTD (full), spatialDWLS, MuSiC, NNLS on the training reference (<= 1,000 cells/type).
# Seeds: one pilot seed (20261101), L'Ecuyer-CMRG; one stream per replicate job, one sub-stream
# per method (ADEMP §8).
#
# Pre-registered checks (written 2026-10-30 before the first run; stop if one fails):
#   C1. truth: rows sum to 1, zero for reference types absent from the spots;
#   C2. every dataset of a job has 500 spots and the same composition (common random numbers);
#   C3. after thinning, no spot exceeds the target (exact thinning since 2026-11-01: the first
#       colon run, with binomial thinning and a margin of 5% + 50 UMI, failed C3 in 9 of 40 jobs
#       by chance; pilot results are not study results);
#   C4. C− μ has the same column sums as E μ;
#   C5. no reference cell comes from the pilot held-out donor.
# Reported: run time per method; failures per method x dataset (ADEMP §7.2); RNG state after each
# method (a method that resets the seed gives the same state in every replicate, ADEMP §8);
# SD(RMSE_i), its upper 95% limit (chi-square, reps - 1 df) and n_sim from that upper limit
# (DECISIONS 2026-10-30) for target MCSE 0.005 / 0.0025 / 0.001; spots at the target depth
# (thinned) and below it, per dataset (exact thinning amendment, DECISIONS 2026-11-02).
# Progress: one line per job start and end in <tissue>_progress.log, next to the summary.
#
# Amendment 2026-11-05 (Luca's email of 2026-11-05, written before the run; rule not yet in DECISIONS):
#   optional 5th argument c = Dirichlet concentration, alpha_k = c * p_k (expected composition still
#   p). Without it alpha_k = K * p_k as before (c = K). Source of c: R6 of the main project, groups of
#   10 neighbouring nuclei (thesis/protocol/R6_dirichlet_concentration_by_roi.csv): colon about 1.
#   With c the outputs go to derived/pilot/<tissue>_c<c>/ and <tissue>_c<c>_* files, so the c = K
#   pilot is untouched; same pilot seed, donors, arms (E, C−, anchor), methods and checks C1–C5.
#   Each replicate file also stores alpha. Expected (pre-registered): the null-baseline RMSE rises to
#   the formula value of 06b_pilot_nsim.R (spots less mixed, so the mean composition guesses worse).
# Amendment 2026-11-06 (DECISIONS 2026-11-06, Luca; before the lymphoid and brain pilots): c = "rule"
#   takes c per tissue and cells level from dirichlet_c() (R/simulate_spots.R; colon 0.68/1.05,
#   lymphoid 1.79/2.05, brain 3.37/21.8 at 2/10 cells); outputs in derived/pilot/<tissue>_crule/.
#   A number still gives one c for both levels, so the colon c = 1 pilot is reproducible.
#
# Usage (working directory = thesis/):  Rscript scripts/06_pilot.R <tissue> [workers] [reps] [mode] [c | rule]
#   "depth" rebuilds the spots from the same seeds and only counts spot depths, without deconvolution.
#   "truth" rebuilds compositions, cells and anchor from the same seeds, checks that the cell truth
#   equals the saved one and adds the UMI-fraction truth (truth_umi) to existing replicate files.

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
workers <- if (length(args) >= 2) as.integer(args[2]) else 8
n_rep <- if (length(args) >= 3) as.integer(args[3]) else 10
mode <- if (length(args) >= 4) args[4] else "full"
stopifnot(mode %in% c("full", "depth", "truth"))
conc <- if (length(args) >= 5) args[5] else NA_character_   # Dirichlet concentration: a number, "rule" or NA (c = K)
stopifnot(is.na(conc) || conc == "rule" || !is.na(suppressWarnings(as.numeric(conc))))
run_id <- if (is.na(conc)) tissue else if (conc == "rule") paste0(tissue, "_crule") else sprintf("%s_c%g", tissue, as.numeric(conc))

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
ref_file <- file.path(data_dir, "derived", "reference", paste0(tissue, "_pilot.rds"))
out_dir <- file.path(data_dir, "derived", "pilot", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(dirname(out_dir), paste0(run_id, "_progress.log"))
depth_file <- file.path(dirname(out_dir), paste0(run_id, "_pilot_depth.csv"))

n_spots <- 500
levels_cells <- c(2, 10)
levels_depth <- c(175, 4900)   # UMI per cell; spot target = n_cells x depth
pilot_seed <- 20261101

R0 <- readRDS(ref_file)
n_train <- lengths(R0$L_by_type[R0$types_spot])          # training cells per spot type
source("R/simulate_spots.R")   # dirichlet_c()
c_for <- function(nc) if (is.na(conc)) length(n_train) else if (conc == "rule") dirichlet_c(tissue, nc) else as.numeric(conc)
alpha_by_n <- setNames(lapply(levels_cells, function(nc) c_for(nc) * n_train / sum(n_train)), levels_cells)   # alpha_k = c * p_k (c = K by default)
if (any(R0$ref_meta$donor_id %in% R0$donors$heldout)) stop("C5: held-out donor in the reference")

grid <- expand.grid(rep = seq_len(n_rep), depth = levels_depth, n_cells = levels_cells)
RNGkind("L'Ecuyer-CMRG")
set.seed(pilot_seed)
s <- .Random.seed
grid$seed <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) { grid$seed[[i]] <- s; s <- parallel::nextRNGStream(s) }

run_job <- function(i) {
  g <- grid[i, ]
  f <- file.path(out_dir, sprintf("c%02d_d%04d_r%02d.rds", g$n_cells, g$depth, g$rep))
  if (mode == "full" && file.exists(f)) return(f)
  job <- sprintf("job %d/%d (cells %d, depth %d, rep %d)", i, nrow(grid), g$n_cells, g$depth, g$rep)
  log_line <- function(msg) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = log_file, append = TRUE)
  log_line(paste("start", job))
  t_job <- proc.time()[["elapsed"]]
  R <- R0
  assign(".Random.seed", g$seed[[1]], envir = .GlobalEnv)
  ts <- R$types_spot
  target <- g$n_cells * g$depth
  comp <- draw_compositions(n_spots, g$n_cells, ts, alpha_by_n[[as.character(g$n_cells)]])
  cells <- draw_cells(comp, R$L_by_type[ts])
  truth <- matrix(0, n_spots, length(R$types_ref), dimnames = list(rownames(comp), R$types_ref))
  truth[, ts] <- comp / g$n_cells
  if (any(abs(rowSums(truth) - 1) > 1e-12)) stop("C1")
  mu_E <- R$mu_E[, ts, drop = FALSE]
  mu_C <- permute_mu(mu_E)
  if (any(abs(colSums(mu_C) - colSums(mu_E)) > 1e-9)) stop("C4")
  held_ok <- R$held_meta$cell_type %in% ts
  anc <- anchor_spots(R$held[, held_ok], R$held_meta$cell_type[held_ok], comp)
  # UMI-fraction truth (secondary estimand, 2026-11-03); draws no random numbers, so the RNG
  # stream of the rest of the job is unchanged. Arms E and C− share it (same cells and L).
  held_used <- R$held[, held_ok][, anc$cells_used, drop = FALSE]
  tu_arm <- umi_truth(cells$spot, cells$type, cells$L, rownames(comp), R$types_ref)
  tu_anc <- umi_truth(anc$spot, R$held_meta$cell_type[held_ok][anc$cells_used],
                      Matrix::colSums(held_used), rownames(comp), R$types_ref)
  truth_umi <- list(E = tu_arm, Cminus = tu_arm, anchor = tu_anc)
  if (any(vapply(truth_umi, function(t) any(abs(rowSums(t) - 1) > 1e-9), logical(1)))) stop("C1 (UMI truth)")
  if (mode == "truth") {   # add the UMI truth to an existing replicate file; estimates untouched
    x <- readRDS(f)
    if (!identical(x$truth, truth)) stop("truth mode: regenerated cell truth differs from the saved one")
    x$truth_umi <- truth_umi
    saveRDS(x, f)
    log_line(paste("truth only, done", job))
    return(f)
  }
  raw_umi <- list()
  thin <- function(d, x) { raw_umi[[d]] <<- Matrix::colSums(x); thin_to_depth(x, target) }   # same RNG order as before
  data <- list(
    E = thin("E", simulate_arm(mu_E, R$phi, cells, rownames(comp))),
    Cminus = thin("Cminus", simulate_arm(mu_C, R$phi, cells, rownames(comp))),
    anchor = thin("anchor", anc$counts))
  for (d in names(data)) {
    if (ncol(data[[d]]) != n_spots) stop("C2: ", d)
    if (max(Matrix::colSums(data[[d]])) > target) stop("C3: ", d)
  }
  depth <- data.frame(tissue = tissue, n_cells = g$n_cells, depth = g$depth, rep = g$rep,
    dataset = names(raw_umi), target = target,
    n_at_target = vapply(raw_umi, function(u) sum(u >= target), numeric(1)),   # thinned to the target
    n_below = vapply(raw_umi, function(u) sum(u < target), numeric(1)),        # left below it
    median_raw_umi = vapply(raw_umi, median, numeric(1)), row.names = NULL)
  if (mode == "depth") { log_line(paste("depth only, done", job)); return(depth) }
  ms <- g$seed[[1]]
  rows <- list(); props <- list()
  for (d in names(data)) {
    x <- data[[d]]
    est <- list()
    for (m in c("NNLS", "RCTD", "MuSiC", "spatialDWLS")) {
      ms <- parallel::nextRNGSubStream(ms)
      assign(".Random.seed", ms, envir = .GlobalEnv)
      r <- tryCatch(switch(m,
        NNLS = run_nnls(x, R$profiles),
        RCTD = run_rctd(x, R$ref, R$ref_meta$cell_type),
        MuSiC = run_music(x, R$ref, R$ref_meta$cell_type, R$ref_meta$donor_id),
        spatialDWLS = run_dwls(x, R$ref, R$ref_meta$cell_type, n_cell = g$n_cells)),
        error = function(e) list(prop = NULL, seconds = NA, error = conditionMessage(e)))
      p <- if (is.null(r$prop)) matrix(NA_real_, n_spots, length(R$types_ref), dimnames = dimnames(truth))
           else complete_prop(r$prop, rownames(truth), R$types_ref)
      est[[m]] <- p
      rows[[length(rows) + 1]] <- data.frame(dataset = d, method = m, seconds = r$seconds,
        error = if (is.null(r$error)) NA else r$error, n_fail = sum(failed_spots(p)),
        rng_after = sum(as.numeric(.Random.seed)), stringsAsFactors = FALSE)
    }
    fs <- lapply(est, failed_spots)
    fs <- fs[!vapply(fs, all, logical(1))]          # a method failing on every spot does not empty the set
    ok <- if (length(fs)) !Reduce(`|`, fs) else rep(FALSE, n_spots)   # spots returned by all methods (ADEMP §7.2)
    for (m in names(est)) {
      k <- which(vapply(rows, function(z) z$dataset == d && z$method == m, logical(1)))
      rows[[k]]$n_common <- sum(ok)
      rows[[k]]$rmse <- if (any(ok)) rmse(est[[m]][ok, ], truth[ok, ]) else NA
      rows[[k]]$jsd <- if (any(ok)) jsd(est[[m]][ok, ], truth[ok, ]) else NA
    }
    props[[d]] <- est
  }
  res <- do.call(rbind, rows)
  res <- cbind(tissue = tissue, n_cells = g$n_cells, depth = g$depth, rep = g$rep, res)
  spot_umi <- vapply(data, function(x) median(Matrix::colSums(x)), numeric(1))
  saveRDS(list(metrics = res, props = props, alpha = alpha_by_n[[as.character(g$n_cells)]], truth = truth, truth_umi = truth_umi, spot_umi = spot_umi, depth = depth,
               anchor_cells_used = length(unique(anc$cells_used)), held_cells = sum(held_ok)), f)
  log_line(sprintf("done %s in %.1f min", job, (proc.time()[["elapsed"]] - t_job) / 60))
  f
}

report_depth <- function(dp) {
  write.csv(dp, depth_file, row.names = FALSE)
  cat("\nSpots at the target depth (thinned) and below it, summed over replicates:\n")
  print(aggregate(cbind(n_at_target, n_below) ~ dataset + depth + n_cells, dp, sum))
  cat("Depth counts written to", depth_file, "\n")
}

cat(sprintf("[%s] pilot (%s): %d jobs (%d reps x %d conditions), %d workers\n",
            tissue, mode, nrow(grid), n_rep, nrow(grid) / n_rep, workers))
cat("Progress log:", log_file, "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("run start, mode %s, %d jobs, %d workers", mode, nrow(grid), workers),
    "\n", file = log_file, append = TRUE)
t0 <- proc.time()[["elapsed"]]
cl <- parallel::makeCluster(workers, setup_strategy = "sequential")   # parallel start-up drops workers on Windows
lib <- .libPaths()
parallel::clusterCall(cl, function(lib) { .libPaths(lib); RNGkind("L'Ecuyer-CMRG"); NULL }, lib)
parallel::clusterExport(cl, c("grid", "R0", "conc", "out_dir", "n_spots", "alpha_by_n", "tissue", "mode", "log_file", "run_job"))
parallel::clusterEvalQ(cl, { for (f in c("R/simulate_spots.R", "R/deconvolve.R", "R/pilot_metrics.R")) source(f); NULL })
done <- parallel::parLapplyLB(cl, seq_len(nrow(grid)), function(i)
  tryCatch(run_job(i), error = function(e) paste("ERROR", i, conditionMessage(e))))
parallel::stopCluster(cl)
errs <- grep("^ERROR", unlist(Filter(is.character, done)), value = TRUE)
if (length(errs)) { cat(errs, sep = "\n"); stop(length(errs), " jobs failed a check") }
cat(sprintf("  done in %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
if (mode == "depth") { report_depth(do.call(rbind, done)); quit(save = "no") }
if (mode == "truth") { cat("UMI truth added to", length(done), "replicate files\n"); quit(save = "no") }

# Summary
out <- lapply(list.files(out_dir, pattern = "^c.*\\.rds$", full.names = TRUE), readRDS)
m <- do.call(rbind, lapply(out, `[[`, "metrics"))
key <- c("n_cells", "depth", "dataset", "method")
agg <- do.call(rbind, lapply(split(m, m[key], drop = TRUE), function(z) data.frame(
  z[1, key], reps = nrow(z), mean_rmse = mean(z$rmse), sd_rmse = sd(z$rmse),
  mean_jsd = mean(z$jsd), sd_jsd = sd(z$jsd), mean_seconds = mean(z$seconds),
  fail_spots = sum(z$n_fail), errors = sum(!is.na(z$error)),
  same_rng_state = length(unique(z$rng_after)) < nrow(z))))
agg$sd_rmse_upper <- agg$sd_rmse * sqrt((agg$reps - 1) / qchisq(0.025, agg$reps - 1))   # upper 95% limit
for (tg in c(0.005, 0.0025, 0.001)) agg[[paste0("n_sim_", tg)]] <- ceiling((agg$sd_rmse_upper / tg)^2)
rownames(agg) <- NULL
write.csv(agg, file.path(dirname(out_dir), paste0(run_id, "_pilot_summary.csv")), row.names = FALSE)
umi <- do.call(rbind, lapply(out, `[[`, "spot_umi"))
cat("\nMedian UMI per spot by dataset:\n"); print(round(apply(umi, 2, median)))
cat(sprintf("Anchor: %d held-out cells of spot types\n", out[[1]]$held_cells))
cat("\nRun time (s) per method, mean over all jobs:\n"); print(round(tapply(m$seconds, m$method, mean, na.rm = TRUE), 1))
cat("\nLargest SD(RMSE_i) per method (over conditions and datasets):\n"); print(signif(cbind(sd = tapply(agg$sd_rmse, agg$method, max), upper95 = tapply(agg$sd_rmse_upper, agg$method, max)), 3))
cat("\nFailed spots and errors per method x dataset:\n"); print(aggregate(cbind(fail_spots, errors) ~ method + dataset, agg, sum))
dp <- do.call(rbind, lapply(out, `[[`, "depth"))   # NULL for jobs run before 2026-11-02
if (!is.null(dp)) report_depth(dp)
if (any(agg$same_rng_state)) cat("\nWARNING: a method leaves the same RNG state in two replicates (seed reset?)\n")
cat("\nSummary written to", file.path(dirname(out_dir), paste0(run_id, "_pilot_summary.csv")), "\n")
