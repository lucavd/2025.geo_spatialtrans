# 06c_pilot_cost.R — cost check of ADEMP §6 from the pilot times.
# Usage: Rscript --vanilla scripts/06c_pilot_cost.R <run_suffix> <n_sim> <workers colon,lymphoid,brain> [days]
#   e.g. Rscript --vanilla scripts/06c_pilot_cost.R crule 169 12,15,5 12
# Amendment 2026-11-17: <run_suffix> and the 5th argument (depth) take one value or one per tissue
# (colon,lymphoid,brain), because the final cost check reads the colon pilot rerun at 175 UMI per cell
# (colon_crule) and the HD pilots of lymphoid and brain (_crule_hd, 62 and 121), with n_sim = 118:
#   Rscript --vanilla scripts/06c_pilot_cost.R crule,crule_hd,crule_hd 118 12,15,5 12 175,62,121
# Optional 5th argument: one depth level (UMI per cell, e.g. 175) to keep, for a design with a single depth
# level (proposal of 2026-11-11: one Visium HD level per tissue; 175 is the pilot level closest to it).
# Reads  $TESI_DATA/derived/pilot/<tissue>_<suffix>_pilot_summary.csv (column mean_seconds, per dataset of 500 spots).
# Assumptions (ADEMP §6): one core per worker (BLAS one thread); the arms not timed in the pilot
# (P, E', O) cost as much as the mean of the timed ones (anchor, E, C-); 6 datasets per condition
# (anchor, E, C-, E', O, P; G would add a 7th, +17%); simulation time and cell2location (GPU) excluded.
# Workers: the three tissues run side by side, so their sum must stay <= 32 cores and their memory
# <= ~85 GB (96 GB cap per user); brain ~16 GB per worker (pilot 10/11).

args <- commandArgs(trailingOnly = TRUE)
suffix  <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else "crule"
n_sim   <- if (length(args) >= 2) as.integer(args[2]) else 169L
workers <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(12L, 15L, 5L)
days    <- if (length(args) >= 4) as.numeric(args[4]) else 12   # weeks 7-8 = 30/11-11/12
depth_only <- if (length(args) >= 5) as.numeric(strsplit(args[5], ",")[[1]]) else NA
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
tissues <- c("colon", "lymphoid", "brain")
names(workers) <- tissues
stopifnot(length(suffix) %in% c(1, 3), length(depth_only) %in% c(1, 3))
suffix <- setNames(rep_len(suffix, 3), tissues); depth_only <- setNames(rep_len(depth_only, 3), tissues)
n_datasets <- 6

rows <- list()
for (t in tissues) {
  f <- file.path(data_dir, "derived", "pilot", sprintf("%s_%s_pilot_summary.csv", t, suffix[[t]]))
  if (!file.exists(f)) { message("missing: ", f); next }
  s <- read.csv(f)
  if (!is.na(depth_only[[t]])) s <- s[s$depth == depth_only[[t]], ]
  if (!nrow(s)) stop("no pilot rows for ", t, " at depth ", depth_only[[t]])
  # mean seconds per method and condition, averaged over the timed datasets
  by_cond <- aggregate(mean_seconds ~ method + n_cells + depth, data = s, FUN = mean)
  per_method <- aggregate(mean_seconds ~ method, data = by_cond, FUN = mean)
  n_cond <- nrow(unique(by_cond[, c("n_cells", "depth")]))
  sec_rep <- sum(by_cond$mean_seconds) * n_datasets   # all conditions, all methods
  cpu_h_rep <- sec_rep / 3600
  cpu_h_tot <- cpu_h_rep * n_sim
  wall_days <- cpu_h_tot / workers[[t]] / 24
  reps_fit  <- floor(days * 24 * workers[[t]] / cpu_h_rep)
  rows[[t]] <- data.frame(tissue = t, run = suffix[[t]], depth = depth_only[[t]], conditions = n_cond,
                          RCTD = per_method$mean_seconds[per_method$method == "RCTD"],
                          MuSiC = per_method$mean_seconds[per_method$method == "MuSiC"],
                          spatialDWLS = per_method$mean_seconds[per_method$method == "spatialDWLS"],
                          NNLS = per_method$mean_seconds[per_method$method == "NNLS"],
                          cpu_h_per_rep = cpu_h_rep, workers = workers[[t]],
                          cpu_h_total = cpu_h_tot, wall_days = wall_days, reps_in_window = reps_fit)
}
out <- do.call(rbind, rows)
cat(sprintf("Cost check: n_sim = %d, %d datasets per condition, window %g days%s\n", n_sim, n_datasets, days,
            if (all(is.na(depth_only))) "" else ", one depth level per tissue (column depth)"))
cat("Seconds per dataset of 500 spots (mean over conditions and timed arms):\n")
print(format(out, digits = 3), row.names = FALSE)
cat(sprintf("\nAll tissues in parallel: %.1f days (slowest tissue); sequential on 32 cores: %.1f days\n",
            max(out$wall_days), sum(out$cpu_h_total) / 32 / 24))
cat("Share of CPU time per method (all tissues):\n")
m <- colSums(out[, c("RCTD", "MuSiC", "spatialDWLS", "NNLS")] * out$cpu_h_total /
             rowSums(out[, c("RCTD", "MuSiC", "spatialDWLS", "NNLS")]))
print(round(m / sum(m), 2))
