
# review/code/recompute.R — ricalcolo indipendente (revisore "code") dai RDS grezzi di S1.3b. Non usa tools/S1.3b_analysis.R.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(parallel))
B <- "/mnt/micron/geo_spatialtrans/S1.3b"; OUTR <- "results/S1.3b/review/code"
L <- function(...) cat(sprintf(...), "\n", sep = "")
nf <- sort(list.files(file.path(B, "null"), pattern = "[.]rds$", full.names = TRUE)); L("null files: %d", length(nf))
rows <- mclapply(nf, function(f) { o <- readRDS(f)
  do.call(rbind, lapply(names(o$res), function(v) { m <- o$res[[v]]; pr <- m$pairs
    gk <- c(V0 = 0, V0g = 0, V1 = 1e-6, V2 = 1e-6, V3 = 1e-7)[[v]]; g2 <- if (gk == 0) 1e-6 else gk
    tau_me <- gk * (o$Preg + o$sumP)
    data.frame(file = basename(f), arch = o$archetype, roi = o$roi_id, model = o$model, rep = o$rep, n = o$n, A = o$A, Preg = o$Preg, sumP = o$sumP,
      holes = o$reg[["n_holes"]], regv = o$reg[["n_vertices"]], v = v, n_out = m$n, ov = m$ov, symd = m$symd, tau_file = m$tau, tau_me = tau_me,
      sum_a = m$sum_a, union_a = m$union_a,
      np = nrow(pr), pov = sum(pr$area > g2 * pr$psum), pov_tauij = if (nrow(pr)) sum(pr$area > pr$tau_ij) else 0L,
      pmax = if (nrow(pr)) max(pr$area) else 0, pairs_ij = if (nrow(pr) && any(pr$area > g2 * pr$psum)) paste(pr$i[pr$area > g2 * pr$psum], pr$j[pr$area > g2 * pr$psum], sep = "/", collapse = ";") else "",
      pairs_area = if (nrow(pr) && any(pr$area > g2 * pr$psum)) paste(signif(pr$area[pr$area > g2 * pr$psum], 6), collapse = ";") else "",
      inv = if (is.null(m$n_invalid)) NA else m$n_invalid, emp = if (is.null(m$n_empty)) NA else m$n_empty,
      coll = if (is.null(m$n_collapsed)) NA else m$n_collapsed,
      types_ok = if (is.null(m$types)) NA else identical(c(m$types), c(o$res$V0$types)),
      cov_n = if (is.null(m$cov)) NA else m$cov[["n_pts"]], cov_anom = if (is.null(m$cov)) NA else m$cov[["n_anom"]],
      cov_far = if (is.null(m$cov)) NA else m$cov[["n_anom_far"]], cov_dmax = if (is.null(m$cov)) NA else m$cov[["dmax"]],
      cov0 = if (is.null(m$cov)) NA else m$cov[["cov0"]], cov2 = if (is.null(m$cov)) NA else m$cov[["cov2"]],
      c5_over = m$c5[["n_over"]], c5_ratio = m$c5[["max_ratio"]], c5_abs = m$c5[["max_abs"]],
      b1_n = m$b1[["n_cells"]], b1_max = m$b1[["max_rel"]], b1_over = m$b1[["n_over"]], flips = m$interior_flips,
      med_area = m$summary[["median_area"]], cv = m$summary[["cv"]], med_eqr = m$summary[["median_eq_r"]], frac_int = m$summary[["frac_interior"]],
      cv_loc = m$summary[["cv_loc"]], ecc = m$summary[["median_ecc_T"]],
      i_diff = if (is.null(m$info)) NA else m$info[["n_pairs_cross_diff"]], i_tj = if (is.null(m$info)) NA else m$info[["n_tjunction"]],
      i_drop = if (is.null(m$info)) NA else m$info[["n_rings_dropped"]],
      t_build = m$time[["build"]], t_union = m$time[["union"]],
      c0_n = o$c0[["c1"]], c0_c1 = o$c0[["c1"]], c0_c2 = o$c0[["c2_overlap"]], c0_sd = o$c0[["c2_symdiff"]])
  })) }, mc.cores = 30)
bad <- vapply(rows, inherits, TRUE, "try-error"); L("read errors: %d", sum(bad))
T <- do.call(rbind, rows); write.csv(T, file.path(OUTR, "null_long.csv"), row.names = FALSE)
L("rows %d; per variant: %s", nrow(T), paste(names(table(T$v)), table(T$v), collapse = " "))
# C-0 against the published S1.3 CSV (independent of old_c123 stored in RDS)
s13 <- read.csv("results/S1.3/S1.3_null_summary.csv"); s13$k <- paste(s13$archetype, s13$roi_id, s13$model, s13$rep)
t0 <- T[T$v == "V0", ]; t0$k <- paste(t0$arch, t0$roi, t0$model, t0$rep); m <- match(t0$k, s13$k)
L("C-0: matched %d; n equal %d; |dc1|<=1e-12 %d; |dc2|<=1e-12 %d; |dsd|<=1e-12 %d; max|dc1| %.3g max|dc2| %.3g",
  sum(!is.na(m)), sum(t0$n == s13$n[m]), sum(abs(t0$c0_c1 - s13$c1[m]) <= 1e-12), sum(abs(t0$c0_c2 - s13$c2_overlap[m]) <= 1e-12),
  sum(abs(t0$c0_sd - s13$c2_symdiff[m]) <= 1e-12), max(abs(t0$c0_c1 - s13$c1[m])), max(abs(t0$c0_c2 - s13$c2_overlap[m])))
f_old <- s13$c2_overlap[m] > 1e-9 | s13$c2_symdiff[m] > 1e-9
f_new <- t0$ov / t0$A > 1e-9 | t0$symd / t0$A > 1e-9
f_new_abs <- abs(t0$ov) / t0$A > 1e-9 | t0$symd / t0$A > 1e-9
L("C-0: S1.3 fails %d; S1.3b V0 fails (signed) %d (abs) %d; identical sets %s / %s", sum(f_old), sum(f_new), sum(f_new_abs), identical(f_old, f_new), identical(f_old, f_new_abs))
L("C-0 fail list: %s", paste(t0$k[f_new_abs], collapse = "; "))
L("V0 ov sign among fails: ov>0 %d ov<0 %d; symd-only %d", sum(t0$ov[f_new_abs] > 0), sum(t0$ov[f_new_abs] < 0), sum(f_new_abs & abs(t0$ov)/t0$A <= 1e-9))
# C-1
for (v in c("V0", "V0g")) { d <- T[T$v == v, ]; L("C-1 %s: nulls |ov| or symd > 1e-9 A: %d", v, sum(abs(d$ov)/d$A > 1e-9 | d$symd/d$A > 1e-9)) }
f0g <- with(T[T$v == "V0g", ], abs(ov)/A > 1e-9 | symd/A > 1e-9); L("C-1 V0g same set as V0: %s", identical(f0g, f_new_abs))
for (v in c("V1", "V2", "V3")) { d <- T[T$v == v, ]
  L("C-1 %s: over tau %d; max|ov|/tau %.3g; max symd/tau %.3g; tau_file==tau_me %s (max rel diff %.2g)", v, sum(abs(d$ov) > d$tau_me | d$symd > d$tau_me),
    max(abs(d$ov)/d$tau_me), max(d$symd/d$tau_me), isTRUE(all.equal(d$tau_file, d$tau_me)), max(abs(d$tau_file/d$tau_me - 1))) }
# C-2
for (v in c("V0", "V0g", "V1", "V2", "V3")) { d <- T[T$v == v, ]; k <- d$pov > 0
  L("C-2 %s: pairs over %d (tau_ij field %d) in %d nulls; nulls: %s", v, sum(d$pov), sum(d$pov_tauij), sum(k), paste(paste(d$arch, d$roi, d$model, d$rep)[k], collapse = "; ")) }
d <- T[T$v == "V0", ]; k <- d$pov > 0; L("C-2 V0 pairs: %s", paste(paste(d$arch, d$roi, d$model, d$rep, d$pairs_ij, d$pairs_area)[k], collapse = " | "))
g <- T[T$v == "V0g", ]; L("C-2 V0g same pairs as V0: %s", identical(g$pairs_ij, d$pairs_ij))
L("C-2 V0 pair nulls also union-failing: %d of %d", sum(k & f_new_abs), sum(k))
L("C-2 V0 total predicate pairs (2********) over all nulls %d; nulls with >=1 predicate pair %d", sum(d$np), sum(d$np > 0))
for (v in c("V1", "V2", "V3")) { e <- T[T$v == v, ]; L("C-2 %s predicate pairs total %d; max pair area %.3g", v, sum(e$np), max(e$pmax)) }
# C-3
for (v in c("V0", "V1", "V2", "V3")) { e <- T[T$v == v, ]
  L("C-3 %s: pts %g; anom %d; far %d; dmax %.3g; cov0 %d; cov2 %d; NA %d", v, sum(e$cov_n), sum(e$cov_anom), sum(e$cov_far), max(e$cov_dmax), sum(e$cov0), sum(e$cov2), sum(is.na(e$cov_n))) }
# C-4
for (v in c("V1", "V2", "V3")) { e <- T[T$v == v, ]
  L("C-4 %s: invalid %d (nulls %d: %s); empty %d; n differs %d; collapsed %d; types changed %d", v, sum(e$inv), sum(e$inv > 0),
    paste(paste(e$arch, e$roi, e$model, e$rep)[e$inv > 0], collapse = "; "), sum(e$emp), sum(e$n_out != e$n), sum(e$coll), sum(!e$types_ok)) }
# C-5
for (v in c("V1", "V2", "V3")) { e <- T[T$v == v, ]; L("C-5 %s: cells over %d; max ratio %.3g; max abs %.3g", v, sum(e$c5_over), max(e$c5_ratio), max(e$c5_abs)) }
# B-1
for (v in c("V1", "V2", "V3")) { e <- T[T$v == v, ]; L("B-1 %s: nulls %d; cells %d; max rel %.3g; cells checked %d", v, sum(e$b1_over > 0), sum(e$b1_over), max(e$b1_max), sum(e$b1_n)) }
# B-2 summaries
k0 <- paste(t0$arch, t0$roi, t0$model, t0$rep)
for (v in c("V1", "V2", "V3")) { e <- T[T$v == v, ]; mm <- match(paste(e$arch, e$roi, e$model, e$rep), k0)
  r <- sapply(c("med_area", "cv", "med_eqr", "frac_int", "cv_loc", "ecc"), function(cc) max(abs(e[[cc]] / t0[[cc]][mm] - 1)))
  L("B-2 %s: %s; interior flips %d", v, paste(names(r), signif(r, 3), collapse = ", "), sum(e$flips)) }
# timings
for (v in c("V1", "V2", "V3")) { e <- T[T$v == v, ]; L("time %s build+union %.1f s", v, sum(e$t_build + e$t_union)) }
# CP-3
fv3 <- T[T$v == "V3", ]; stopifnot(identical(paste(fv3$arch, fv3$roi, fv3$model, fv3$rep), k0))
for (cc in list(c("n", "V0"), c("holes", "V0"), c("regv", "V0"), c("i_diff", "V3"), c("i_tj", "V3"))) { x <- T[T$v == cc[2], cc[1]]
  L("CP-3 %s: median fail %g pass %g p %.2g", cc[1], median(x[f_new_abs]), median(x[!f_new_abs]), wilcox.test(x ~ f_new_abs)$p.value) }
L("CP-3 V3 n_pairs_cross_diff among fails: %s", paste(fv3$i_diff[f_new_abs], collapse = ","))
L("CP-3 V3 i_diff == 0 among pass: %d of %d; >=3 among pass %d", sum(fv3$i_diff[!f_new_abs] == 0), sum(!f_new_abs), sum(fv3$i_diff[!f_new_abs] >= 3))
