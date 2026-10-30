# 06_pilot.R — declared pilot (PROTOCOLLO §6 week 2; ADEMP §10). Its numbers never enter the
# thesis results: it only fills the [PILOT] placeholders (run time, failures, SD(RMSE_i) for n_sim,
# seed check, reuse of held-out cells).
#
# Design (ADEMP §10, proposal 2026-10-29; levels Pierangela 2026-10-30):
#   pilot donors from 03_reference.R <tissue> pilot (drawn inside training, seed 20261030);
#   arms E and C− (one new gene permutation per replicate), plus the real anchor (built from the
#   pilot held-out donor), because n_sim takes the largest SD across the design (ADEMP §6) and
#   the anchor is part of every condition;
#   cells per spot 5 / 15; depth 10,000 / 20,000 UMI per spot (spots above are thinned);
#   composition Dirichlet(alpha = 1) over the spot types, cells drawn multinomially;
#   500 spots per dataset; 10 replicates per condition;
#   methods RCTD (full), spatialDWLS, MuSiC, NNLS on the training reference (<= 1,000 cells/type).
# Seeds: one pilot seed (20261101), L'Ecuyer-CMRG; one stream per replicate job, one sub-stream
# per method (ADEMP §8).
#
# Pre-registered checks (written 2026-10-30 before the first run; stop if one fails):
#   C1. truth: rows sum to 1, zero for reference types absent from the spots;
#   C2. every dataset of a job has 500 spots and the same composition (common random numbers);
#   C3. after thinning, no spot exceeds the target by more than 5% (binomial thinning is not exact);
#   C4. C− μ has the same column sums as E μ;
#   C5. no reference cell comes from the pilot held-out donor.
# Reported: run time per method; failures per method x dataset (ADEMP §7.2); RNG state after each
# method (a method that resets the seed gives the same state in every replicate, ADEMP §8);
# SD(RMSE_i), its upper 95% limit (chi-square, reps - 1 df) and n_sim from that upper limit
# (DECISIONS 2026-10-30) for target MCSE 0.005 / 0.0025 / 0.001.
#
# Usage (working directory = thesis/):  Rscript scripts/06_pilot.R <tissue> [workers] [reps]

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
workers <- if (length(args) >= 2) as.integer(args[2]) else 8
n_rep <- if (length(args) >= 3) as.integer(args[3]) else 10

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
ref_file <- file.path(data_dir, "derived", "reference", paste0(tissue, "_pilot.rds"))
out_dir <- file.path(data_dir, "derived", "pilot", tissue)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_spots <- 500
alpha <- 1
levels_cells <- c(5, 15)
levels_depth <- c(10000, 20000)
pilot_seed <- 20261101

R0 <- readRDS(ref_file)
if (any(R0$ref_meta$donor_id %in% R0$donors$heldout)) stop("C5: held-out donor in the reference")

grid <- expand.grid(rep = seq_len(n_rep), depth = levels_depth, n_cells = levels_cells)
RNGkind("L'Ecuyer-CMRG")
set.seed(pilot_seed)
s <- .Random.seed
grid$seed <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) { grid$seed[[i]] <- s; s <- parallel::nextRNGStream(s) }

run_job <- function(i) {
  g <- grid[i, ]
  f <- file.path(out_dir, sprintf("c%02d_d%05d_r%02d.rds", g$n_cells, g$depth, g$rep))
  if (file.exists(f)) return(f)
  R <- R0
  assign(".Random.seed", g$seed[[1]], envir = .GlobalEnv)
  ts <- R$types_spot
  comp <- draw_compositions(n_spots, g$n_cells, ts, alpha)
  cells <- draw_cells(comp, R$L_by_type[ts])
  truth <- matrix(0, n_spots, length(R$types_ref), dimnames = list(rownames(comp), R$types_ref))
  truth[, ts] <- comp / g$n_cells
  if (any(abs(rowSums(truth) - 1) > 1e-12)) stop("C1")
  mu_E <- R$mu_E[, ts, drop = FALSE]
  mu_C <- permute_mu(mu_E)
  if (any(abs(colSums(mu_C) - colSums(mu_E)) > 1e-9)) stop("C4")
  held_ok <- R$held_meta$cell_type %in% ts
  anc <- anchor_spots(R$held[, held_ok], R$held_meta$cell_type[held_ok], comp)
  data <- list(
    E = thin_to_depth(simulate_arm(mu_E, R$phi, cells, rownames(comp)), g$depth),
    Cminus = thin_to_depth(simulate_arm(mu_C, R$phi, cells, rownames(comp)), g$depth),
    anchor = thin_to_depth(anc$counts, g$depth))
  for (d in names(data)) {
    if (ncol(data[[d]]) != n_spots) stop("C2: ", d)
    if (max(Matrix::colSums(data[[d]])) > 1.05 * g$depth) stop("C3: ", d)
  }
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
    ok <- !Reduce(`|`, lapply(est, failed_spots))   # spots returned by all methods (ADEMP §7.2)
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
  saveRDS(list(metrics = res, spot_umi = spot_umi, anchor_cells_used = length(unique(anc$cells_used)),
               held_cells = sum(held_ok)), f)
  f
}

cat(sprintf("[%s] pilot: %d jobs (%d reps x %d conditions), %d workers\n",
            tissue, nrow(grid), n_rep, nrow(grid) / n_rep, workers))
t0 <- proc.time()[["elapsed"]]
cl <- parallel::makeCluster(workers)
lib <- .libPaths()
parallel::clusterCall(cl, function(lib) { .libPaths(lib); RNGkind("L'Ecuyer-CMRG"); NULL }, lib)
parallel::clusterExport(cl, c("grid", "R0", "out_dir", "n_spots", "alpha", "tissue"))
parallel::clusterEvalQ(cl, { for (f in c("R/simulate_spots.R", "R/deconvolve.R", "R/pilot_metrics.R")) source(f); NULL })
done <- parallel::parLapplyLB(cl, seq_len(nrow(grid)), function(i)
  tryCatch(run_job(i), error = function(e) paste("ERROR", i, conditionMessage(e))))
parallel::stopCluster(cl)
errs <- grep("^ERROR", unlist(done), value = TRUE)
if (length(errs)) { cat(errs, sep = "\n"); stop(length(errs), " jobs failed a check") }
cat(sprintf("  done in %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))

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
write.csv(agg, file.path(dirname(out_dir), paste0(tissue, "_pilot_summary.csv")), row.names = FALSE)
umi <- do.call(rbind, lapply(out, `[[`, "spot_umi"))
cat("\nMedian UMI per spot by dataset:\n"); print(round(apply(umi, 2, median)))
cat(sprintf("Anchor: %d held-out cells of spot types\n", out[[1]]$held_cells))
cat("\nRun time (s) per method, mean over all jobs:\n"); print(round(tapply(m$seconds, m$method, mean, na.rm = TRUE), 1))
cat("\nLargest SD(RMSE_i) per method (over conditions and datasets):\n"); print(signif(cbind(sd = tapply(agg$sd_rmse, agg$method, max), upper95 = tapply(agg$sd_rmse_upper, agg$method, max)), 3))
cat("\nFailed spots and errors per method x dataset:\n"); print(aggregate(cbind(fail_spots, errors) ~ method + dataset, agg, sum))
if (any(agg$same_rng_state)) cat("\nWARNING: a method leaves the same RNG state in two replicates (seed reset?)\n")
cat("\nSummary written to", file.path(dirname(out_dir), paste0(tissue, "_pilot_summary.csv")), "\n")
