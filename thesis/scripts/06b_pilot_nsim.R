# 06b_pilot_nsim.R — pilot tables for n_sim (ADEMP §6, §7.2; DECISIONS 2026-11-03). Pilot numbers
# never enter the thesis results: they fix n_sim and document the failure rule.
#
# Input: the replicate files written by 06_pilot.R (<TESI_DATA>/derived/pilot/<tissue>/c*_d*_r*.rds)
# and the pilot reference (<TESI_DATA>/derived/reference/<tissue>_pilot.rds).
# Rule of the main analysis: failed spots get the null-baseline proportions (share of each type
# among all training cells, all reference types); the null baseline itself is reported as a
# reference value, scored on every spot.
# Measures per replicate: RMSE of each method and of the null baseline; failure rate; paired
# difference D = RMSE(m1) − RMSE(m2) within the replicate. Across replicates: mean, SD and the
# upper 95% limit of the SD (chi-square, reps − 1 df; DECISIONS 2026-10-30).
# n_sim = (4 × SD_upper / δ)², from 2 × MCSE(D) ≤ δ/2, with SD_upper the largest over pairs,
# conditions and datasets (DECISIONS 2026-11-03).
# RMSE relative to the null baseline (DECISIONS 2026-11-04, Luca; ADEMP §5.2; pre-registered
# 2026-11-05 before the run): rel_null = mean RMSE of the method / mean RMSE of the null baseline,
# same condition x dataset (ratio of means); rel_null_umi the same on the UMI-share truth; and the
# paired difference RMSE_i(method) − RMSE_i(null) within replicate, mean and MCSE = SD/√reps, so a
# ratio below 1 is read together with its Monte Carlo error. Check (second route): equal (< 1e-9)
# to the one-off computation of 2026-11-04 (colon_pilot_rmse_relative_null.csv) when that file
# exists. Correction 2026-11-05 after the first run: the check was written against its column
# `ratio`, which is the mean of the per-replicate ratios, not the ratio of means of ADEMP §5.2
# (largest difference 0.0016; its `ratio_umi` is already a ratio of means); it now compares rmse / null
# of that file, which is the ratio of means, and prints the largest difference from the mean of ratios.
#
# Output (next to the pilot summary, outside git):
#   <tissue>_pilot_rmse_methods.csv   one row per condition x dataset x method (+ null baseline)
#   <tissue>_pilot_rmse_pairs.csv     one row per condition x dataset x pair of methods
#   <tissue>_pilot_nsim.csv           n_sim and share of pairs with |mean D| >= δ, per δ
#   <tissue>_pilot_tables.md          the same, formatted (mean (SD)); Tables 5–7 added 2026-11-03
#   <tissue>_pilot_type_diagnostics.csv  per replicate: bias per spot type, proportion on absent types
#
# Usage (working directory = thesis/):  Rscript scripts/06b_pilot_nsim.R <tissue>

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
source("R/pilot_metrics.R")

pilot_dir <- file.path(data_dir, "derived", "pilot")
files <- list.files(file.path(pilot_dir, tissue), pattern = "^c\\d+_d\\d+_r\\d+\\.rds$", full.names = TRUE)
stopifnot("no pilot replicate files" = length(files) > 0)
R0 <- readRDS(file.path(data_dir, "derived", "reference", paste0(tissue, "_pilot.rds")))
p_null <- null_props(lengths(R0$L_by_type))
stopifnot(abs(sum(p_null) - 1) < 1e-12, setequal(names(p_null), R0$types_ref))
# UMI-fraction estimand (secondary, 2026-11-03): null baseline = share of all training UMI per type
p_null_umi <- vapply(R0$L_by_type, sum, numeric(1)); p_null_umi <- p_null_umi / sum(p_null_umi)
methods <- c("RCTD", "spatialDWLS", "MuSiC", "NNLS")
deltas <- c(0.005, 0.01, 0.02, 0.03)
up <- function(n) sqrt((n - 1) / qchisq(0.025, n - 1))   # upper 95% limit factor for an SD

# Per replicate
rows <- list(); brows <- list()
ts <- R0$types_spot; absent <- setdiff(R0$types_ref, ts)
for (f in files) {
  x <- readRDS(f); m0 <- x$metrics[1, ]
  for (ds in names(x$props)) {
    tr <- x$truth
    tu <- if (is.null(x$truth_umi)) NULL else x$truth_umi[[ds]][, colnames(tr)]
    for (m in names(x$props[[ds]])) {
      est <- x$props[[ds]][[m]][, colnames(tr)]
      bad <- failed_spots(est); filled <- fill_failed(est, p_null)
      stopifnot(all(abs(rowSums(filled) - 1) < 1e-6), identical(filled[!bad, ], est[!bad, ]))
      rows[[length(rows) + 1]] <- data.frame(n_cells = m0$n_cells, depth = m0$depth, rep = m0$rep,
        dataset = ds, method = m, rmse = rmse(filled, tr), fail_rate = mean(bad),
        rmse_umi = if (is.null(tu)) NA else rmse(fill_failed(est, p_null_umi), tu))
      # Diagnostics (2026-11-03): mean bias per spot type (cell truth, returned spots) and mean
      # proportion given to reference types absent from the spots
      if (any(!bad)) brows[[length(brows) + 1]] <- data.frame(n_cells = m0$n_cells, depth = m0$depth,
        rep = m0$rep, dataset = ds, method = m, type = c(ts, absent),
        value = c(colMeans(est[!bad, ts, drop = FALSE] - tr[!bad, ts, drop = FALSE]),
                  colMeans(est[!bad, absent, drop = FALSE])),
        kind = rep(c("bias", "absent"), c(length(ts), length(absent))))
    }
    nb <- matrix(p_null[colnames(tr)], nrow(tr), ncol(tr), byrow = TRUE)
    rows[[length(rows) + 1]] <- data.frame(n_cells = m0$n_cells, depth = m0$depth, rep = m0$rep,
      dataset = ds, method = "null baseline", rmse = rmse(nb, tr), fail_rate = 0,
      rmse_umi = if (is.null(tu)) NA else rmse(matrix(p_null_umi[colnames(tr)], nrow(tr), ncol(tr), byrow = TRUE), tu))
  }
}
d <- do.call(rbind, rows)
key <- c("n_cells", "depth", "dataset")

# Single methods
sm <- do.call(rbind, lapply(split(d, d[c(key, "method")], drop = TRUE), function(z) data.frame(
  z[1, c(key, "method")], reps = nrow(z), mean_rmse = mean(z$rmse), sd_rmse = sd(z$rmse),
  sd_rmse_upper = sd(z$rmse) * up(nrow(z)), mean_fail_rate = mean(z$fail_rate), sd_fail_rate = sd(z$fail_rate),
  mean_rmse_umi = mean(z$rmse_umi), sd_rmse_umi = sd(z$rmse_umi))))

# RMSE relative to the null baseline (DECISIONS 2026-11-04): ratio of means, and paired difference
# method − null within replicate with its MCSE
nbl <- d[d$method == "null baseline", ]
rk <- function(z) paste(z$n_cells, z$depth, z$dataset, z$rep)
dn <- d$rmse - nbl$rmse[match(rk(d), rk(nbl))]
dd <- do.call(rbind, lapply(split(data.frame(d[c(key, "method")], dn), d[c(key, "method")], drop = TRUE),
  function(z) data.frame(z[1, c(key, "method")], diff_null = mean(z$dn), mcse_diff_null = sd(z$dn) / sqrt(nrow(z)))))
gk <- function(z) paste(z$n_cells, z$depth, z$dataset)
nb_sm <- sm[sm$method == "null baseline", ]
sm$rel_null <- sm$mean_rmse / nb_sm$mean_rmse[match(gk(sm), gk(nb_sm))]
sm$rel_null_umi <- sm$mean_rmse_umi / nb_sm$mean_rmse_umi[match(gk(sm), gk(nb_sm))]
j <- match(paste(gk(sm), sm$method), paste(gk(dd), dd$method))
sm$diff_null <- dd$diff_null[j]; sm$mcse_diff_null <- dd$mcse_diff_null[j]
old <- file.path(pilot_dir, paste0(tissue, "_pilot_rmse_relative_null.csv"))
if (file.exists(old)) {
  o <- read.csv(old); i <- match(paste(o$n_cells, o$depth, o$dataset, o$method), paste(gk(sm), sm$method))
  stopifnot("rel_null differs from the 2026-11-04 computation" = !anyNA(i),
            max(abs(sm$rel_null[i] - o$rmse / o$null)) < 1e-9)
  cat(sprintf("Relative RMSE: equal to rmse / null of %s (second route); largest difference from its mean of per-replicate ratios: cells %.4f, UMI %.4f\n",
              basename(old), max(abs(sm$rel_null[i] - o$ratio)), max(abs(sm$rel_null_umi[i] - o$ratio_umi))))
}

# Paired differences (methods only, same replicate)
dm <- d[d$method %in% methods, ]
pairs <- combn(methods, 2, simplify = FALSE)
pm <- do.call(rbind, lapply(split(dm, dm[key], drop = TRUE), function(z) do.call(rbind, lapply(pairs, function(p) {
  a <- z[z$method == p[1], ]; b <- z[z$method == p[2], ]
  b <- b[match(a$rep, b$rep), ]; D <- a$rmse - b$rmse
  data.frame(z[1, key], pair = paste(p[1], "-", p[2]), reps = length(D), mean_diff = mean(D),
    sd_diff = sd(D), sd_diff_upper = sd(D) * up(length(D)))
}))))

# n_sim per δ; share of pairs above δ in the arms of the main analysis present in the pilot (E, anchor)
sd_max <- max(pm$sd_diff_upper)
main <- pm$dataset %in% c("E", "anchor")
ns <- data.frame(delta = deltas, sd_diff_upper_max = sd_max, n_sim = ceiling((4 * sd_max / deltas)^2),
  share_pairs_ge_delta = sapply(deltas, function(dl) mean(abs(pm$mean_diff[main]) >= dl)),
  n_pairs = sum(main))

ord <- function(z) z[order(z$n_cells, z$depth, match(z$dataset, c("E", "Cminus", "anchor"))), ]
sm <- ord(sm); pm <- ord(pm); rownames(sm) <- rownames(pm) <- NULL
out <- function(suffix) file.path(pilot_dir, paste0(tissue, "_pilot_", suffix))
write.csv(sm, out("rmse_methods.csv"), row.names = FALSE)
write.csv(pm, out("rmse_pairs.csv"), row.names = FALSE)
write.csv(ns, out("nsim.csv"), row.names = FALSE)

# Formatted tables: rows condition x method (or pair), columns E / C− / anchor, cells "mean (SD)"
fmt <- function(m, s) sprintf("%.3f (%.3f)", m, s)
md_table <- function(z, rowvar, mcol, scol, rowlabels) {
  cond <- unique(z[c("n_cells", "depth")]); lines <- character()
  for (i in seq_len(nrow(cond))) for (r in rowlabels) {
    cell <- function(ds) { w <- z[z$n_cells == cond$n_cells[i] & z$depth == cond$depth[i] & z$dataset == ds & z[[rowvar]] == r, ]
      if (!nrow(w)) "" else if (is.null(scol)) sprintf("%.2f", w[[mcol]]) else fmt(w[[mcol]], w[[scol]]) }
    lines <- c(lines, sprintf("| %d cells, %s UMI | %s | %s | %s | %s |", cond$n_cells[i],
      format(cond$depth[i], big.mark = ","), r, cell("E"), cell("Cminus"), cell("anchor")))
  }
  c("| Condition | %s | E | C− | Anchor |", "|---|---|---|---|---|", lines)
}
t1 <- md_table(sm, "method", "mean_rmse", "sd_rmse", c(methods, "null baseline")); t1[1] <- sprintf(t1[1], "Method")
t2 <- md_table(sm[sm$method %in% methods, ], "method", "mean_fail_rate", "sd_fail_rate", methods); t2[1] <- sprintf(t2[1], "Method")
t3 <- md_table(pm, "pair", "mean_diff", "sd_diff", unique(pm$pair)); t3[1] <- sprintf(t3[1], "Pair (m1 − m2)")
t4 <- c("| δ | n_sim | Pairs with mean difference ≥ δ (E and anchor) |", "|---|---|---|",
  sprintf("| %.3f | %d | %.0f%% of %d |", ns$delta, ns$n_sim, 100 * ns$share_pairs_ge_delta, ns$n_pairs))
t5 <- md_table(sm, "method", "mean_rmse_umi", "sd_rmse_umi", c(methods, "null baseline")); t5[1] <- sprintf(t5[1], "Method")
# Tables 6–7: mean over replicates, one row per condition x dataset (E, anchor) x method, one column per type
bd <- do.call(rbind, brows); bd <- bd[bd$dataset %in% c("E", "anchor"), ]
type_table <- function(kind, types) {
  z <- aggregate(value ~ n_cells + depth + dataset + method + type, bd[bd$kind == kind, ], mean)
  w <- reshape(z, idvar = c("n_cells", "depth", "dataset", "method"), timevar = "type", direction = "wide")
  names(w) <- sub("^value[.]", "", names(w))
  if (kind == "absent") w$total <- rowSums(w[types])
  w <- w[order(w$n_cells, w$depth, match(w$dataset, c("E", "anchor")), match(w$method, methods)), ]
  cols <- c(types, if (kind == "absent") "total")
  c(paste0("| Condition | Dataset | Method | ", paste(cols, collapse = " | "), " |"),
    paste0("|---|---|---|", strrep("---|", length(cols))),
    sprintf("| %d cells, %s UMI | %s | %s | %s |", w$n_cells, format(w$depth, big.mark = ","), w$dataset, w$method,
      apply(w[cols], 1, function(v) paste(sprintf("%.3f", v), collapse = " | "))))
}
t6 <- type_table("bias", ts); t7 <- type_table("absent", absent)
t8 <- md_table(sm[sm$method %in% methods, ], "method", "rel_null", NULL, methods); t8[1] <- sprintf(t8[1], "Method")
t8b <- md_table(sm[sm$method %in% methods, ], "method", "diff_null", "mcse_diff_null", methods); t8b[1] <- sprintf(t8b[1], "Method")
t8c <- md_table(sm[sm$method %in% methods, ], "method", "rel_null_umi", NULL, methods); t8c[1] <- sprintf(t8c[1], "Method")
write.csv(bd, out("type_diagnostics.csv"), row.names = FALSE)
reps <- unique(sm$reps); n_spots <- nrow(readRDS(files[1])$truth)
md <- c(sprintf("# Pilot %s: RMSE, failures and n_sim", tissue), "",
  sprintf("Declared pilot (ADEMP §10): %d replicates x %d spots per condition and dataset; not study results. Generated by `scripts/06b_pilot_nsim.R` on %s.",
    reps, n_spots, format(Sys.Date())), "",
  "Failed spots get the null-baseline proportions (share of each type among all training cells, all reference types; DECISIONS 2026-11-03). The null baseline row scores these proportions on every spot.", "",
  "## Table 1. RMSE per method, mean (SD) over replicates", "", t1, "",
  "## Table 2. Failure rate per method (fraction of spots), mean (SD) over replicates", "", t2, "",
  "## Table 3. Paired difference in RMSE within replicate (m1 − m2), mean (SD) over replicates", "",
  "A negative value means that m1 has the lower error.", "", t3, "",
  sprintf("## Table 4. Number of replicates: n_sim = (4 x SD_upper / δ)², SD_upper = %.4f (largest upper 95%% limit of SD(D) over pairs, conditions and datasets)", sd_max), "", t4, "",
  "## Table 5. RMSE against the UMI-fraction truth (secondary estimand, 2026-11-03), mean (SD) over replicates", "",
  "Truth = share of the spot's UMI from each type (library sizes of its cells); failed spots and the null baseline use the share of all training UMI per type. RCTD weights are defined on this scale (Cable et al. 2022).", "",
  if (all(is.na(sm$mean_rmse_umi))) "Not available: run `06_pilot.R <tissue> <workers> <reps> truth` first." else t5, "",
  "## Table 6. Mean bias per spot type (estimate − cell truth), returned spots, mean over replicates (E and anchor)", "",
  "Positive = the method gives the type more than its true share of cells.", "", t6, "",
  "## Table 7. Mean proportion given to reference types absent from the spots (true share 0), returned spots, mean over replicates (E and anchor)", "", t7, "",
  "## Table 8. RMSE relative to the null baseline (DECISIONS 2026-11-04): mean RMSE of the method / mean RMSE of the null baseline, same condition and dataset", "",
  "Below 1 the method adds information beyond the expected composition.", "", t8, "",
  "## Table 8b. Paired difference RMSE(method) − RMSE(null baseline) within replicate, mean (MCSE)", "", t8b, "",
  "## Table 8c. As Table 8, on the UMI-share truth (secondary estimand)", "",
  if (all(is.na(sm$rel_null_umi))) "Not available." else t8c, "")
writeLines(md, out("tables.md"), useBytes = TRUE)

cat(sprintf("Largest SD_upper: paired difference %.4f, single method %.4f\n", sd_max,
  max(sm$sd_rmse_upper[sm$method %in% methods])))
print(ns, row.names = FALSE)
cat("Written:", paste(basename(c(out("rmse_methods.csv"), out("rmse_pairs.csv"), out("nsim.csv"), out("tables.md")))), "\n")
