
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(parallel))
B <- "/mnt/micron/geo_spatialtrans/S1.3b"; OUTR <- "results/S1.3b/review/code"
L <- function(...) cat(sprintf(...), "\n", sep = "")
T <- read.csv(file.path(OUTR, "null_long.csv")); t0 <- T[T$v == "V0", ]; k0 <- paste(t0$arch, t0$roi, t0$model, t0$rep)
fail <- abs(t0$ov)/t0$A > 1e-9 | t0$symd/t0$A > 1e-9; pairnull <- t0$pov > 0
# ---- C-6
gf <- list.files(file.path(B, "grid"), pattern = "[.]rds$", full.names = TRUE); L("grid files %d", length(gf))
G <- lapply(gf, readRDS); gk <- vapply(G, function(o) paste(o$archetype, o$roi_id, o$model, o$rep), "")
gs <- vapply(G, function(o) o$set, "")
L("grid sets: fail %d pass %d; fail set == S1.3 fails: %s", sum(gs == "fail"), sum(gs == "pass"), setequal(gk[gs == "fail"], k0[fail]))
L("grid pass per archetype: %s", paste(names(table(substr(gk[gs == "pass"], 1, 2))), table(substr(gk[gs == "pass"], 1, 2)), collapse = " "))
L("grid pass set overlaps S1.3 fails: %d", sum(gk[gs == "pass"] %in% k0[fail]))
L("pair-error nulls (C-2 V0) covered by C-6 grid: %d of %d (%s)", sum(k0[pairnull] %in% gk), sum(pairnull), paste(k0[pairnull][k0[pairnull] %in% gk], collapse = ";"))
for (v in c("V0", "V1", "V2", "V3")) { M <- do.call(rbind, lapply(G, function(o) o[[v]]))
  L("C-6 %s: cell_viol %d bin_viol %d; max cell %.3g; max bin %.3g; pieces %d; n_bnd %d; time %.0f s", v, sum(M[, "cell_viol"]), sum(M[, "bin_viol"]),
    max(M[, "cell_max"]), max(M[, "bin_max"]), sum(M[, "n_pieces"]), sum(M[, "n_bnd"]), sum(M[, "elapsed_s"])) }
M0 <- do.call(rbind, lapply(G, function(o) o$V0)); M2 <- do.call(rbind, lapply(G, function(o) o$V2))
L("C-6 V0 vs V2 pieces identical counts: %d/30", sum(M0[, "n_pieces"] == M2[, "n_pieces"]))
# ---- mutants
mu <- readRDS(file.path(B, "mutants", "mutants.rds")); print(mu)
# ---- CP-1
cp <- readRDS(file.path(B, "cp1", "cp1.rds")); k1 <- cp$case1; k2 <- cp$case2
L("CP-1: a_v0 %.6f relate %s a_cp1a %.3g clipped %s types %s", k1$a_v0, k1$relate, k1$a_cp1a, paste(k1$clipped, collapse = ","), paste(k1$types, collapse = ","))
print(k1$d0_info); L("d0_d: %s; d0_dtj: %s; diff_dist_to_region: %s", paste(k1$d0_d, collapse = ","), paste(k1$d0_dtj, collapse = ","), paste(k1$diff_dist_to_region, collapse = ","))
L("placebo: n_far1 %d n_far2 %d; areas range %.6f..%.6f; |area - a_v0| <= 0.01: %d/20; valid all %s; moved %s", k1$n_far1, k1$n_far2, min(k1$placebo$area), max(k1$placebo$area),
  sum(abs(k1$placebo$area - k1$a_v0) <= 0.01), all(k1$placebo$valid), paste(unique(k1$placebo$moved), collapse = ","))
L("placebo areas identical to a_v0 bitwise: %d/20; max |diff| %.3g", sum(k1$placebo$area == k1$a_v0), max(abs(k1$placebo$area - k1$a_v0)))
L("CP-1c: union_deficit_v0 %.4f symd_v0 %s; cp1c_rel_ov %.3g cp1c_rel_symd %s (len %d)", k2$union_deficit_v0, paste(k2$symd_v0, collapse = ","), k2$cp1c_rel_ov,
  paste(k2$cp1c_rel_symd, collapse = ","), length(k2$cp1c_rel_symd)); print(k2$d0_info); print(k2$v3_info); L("diff_on_region_boundary: %s", paste(k2$diff_on_region_boundary, collapse = ","))
# ---- CP-2
cf <- list.files(file.path(B, "cp2"), pattern = "[.]rds$", full.names = TRUE); L("cp2 files %d", length(cf))
C2 <- do.call(rbind, lapply(cf, function(f) { o <- readRDS(f); r <- o$rows; if (is.null(r$error)) r$error <- NA_character_
  cbind(key = paste(o$archetype, o$roi_id, o$model, o$rep), A = o$A, sumP = o$sumP, Preg = o$Preg, r) }))
sel <- read.csv(file.path(B, "cp2", "selection.csv")); sk <- paste(sel$archetype, sel$roi_id, sel$model, sel$rep)
L("cp2 selection %d, files match %s, in S1.3 fail set %d, per arch %s, overlap pair-error nulls %d", length(sk), setequal(sk, unique(C2$key)), sum(sk %in% k0[fail]),
  paste(table(substr(sk, 1, 2)), collapse = ","), sum(sk %in% k0[pairnull]))
C2$exc <- !is.na(C2$error)
thr0 <- 1e-9 * C2$A + sqrt(2) * C2$delta * C2$sumP
C2$silent <- !C2$exc & ifelse(C2$variant == "V0", abs(C2$ov) > thr0 | C2$symd > thr0, abs(C2$ov) > C2$tau | C2$symd > C2$tau)
C2$fail_posthoc <- C2$exc | C2$silent
for (v in c("V0", "V1", "V2", "V3")) { d <- C2[C2$variant == v, ]
  L("CP-2 %s: per delta (1e-12/1e-10/1e-8) n %s; exceptions %s; silent fails %s; fail rate post hoc %s; silent rate %s; pairs_over silent %s",
    v, paste(tapply(d$exc, d$delta, length), collapse = "/"), paste(tapply(d$exc, d$delta, sum), collapse = "/"), paste(tapply(d$silent, d$delta, sum), collapse = "/"),
    paste(sprintf("%.3f", tapply(d$fail_posthoc, d$delta, mean)), collapse = "/"), paste(sprintf("%.3f", tapply(d$silent, d$delta, mean)), collapse = "/"),
    paste(tapply(d$n_pairs_over, d$delta, function(x) sum(x, na.rm = TRUE)), collapse = "/")) }
d <- C2[C2$variant == "V0", ]
L("CP-2 V0 exceptions with invalid input: %d of %d; invalid inputs total %s", sum(d$exc & d$n_invalid_input > 0), sum(d$exc), paste(tapply(d$n_invalid_input > 0, d$delta, sum), collapse = "/"))
L("CP-2 V0 max |ov|/thr and symd/thr among non-exception rows: %.3g %.3g", max(abs(d$ov[!d$exc]) / thr0[C2$variant == "V0"][!d$exc]), max(d$symd[!d$exc] / thr0[C2$variant == "V0"][!d$exc]))
L("CP-2 V0 max |ov|/(1e-9 A) among non-exception rows (S1.3 threshold only): %.3g", max(abs(d$ov[!d$exc]) / (1e-9 * d$A[!d$exc])))
L("CP-2 V0 error messages: %s", paste(unique(substr(d$error[d$exc], 1, 90)), collapse = " || "))
for (v in c("V1", "V2", "V3")) { e <- C2[C2$variant == v, ]; L("CP-2 %s max |ov|/tau %.3g max symd/tau %.3g", v, max(abs(e$ov)/e$tau), max(e$symd/e$tau)) }
nrm <- 15/1800; r_post <- tapply(d$fail_posthoc, d$delta, mean); r_sil <- tapply(d$silent, d$delta, mean)
L("CP-2 V0 verdict with post hoc rule: %s; without (exceptions excluded or not FAIL): %s", any(r_post > nrm), any(r_sil > nrm))
L("binomial P(X>=1 | n=60, p=15/1800) = %.3f; P(X>=1 | n=180) = %.3f", 1 - dbinom(0, 60, nrm), 1 - dbinom(0, 180, nrm))
# ---- B-2: verdict files
ref <- read.csv("results/S1.3/S1.3_verdicts.csv"); L("S1.3 verdict rows %d; B rows %d; B predictions (confirmed non NA) %d", nrow(ref), sum(ref$section == "B"), sum(ref$section == "B" & !is.na(ref$confirmed)))
for (v in c("V0", "V1", "V2", "V3")) { x <- read.csv(file.path(B, "b2", sprintf("verdicts_%s.csv", v)))
  same_n <- nrow(x) == nrow(ref); same_key <- same_n && identical(paste(x$section, x$id, x$archetype, x$model), paste(ref$section, ref$id, ref$archetype, ref$model))
  kb <- ref$section == "B"; kp <- kb & !is.na(ref$confirmed)
  L("B-2 %s: rows %d keys same %s; all cols identical %s; B outcome same %d/%d; B pred confirmed same %d/%d; B value same %d/%d; non-B rows identical %d/%d", v, nrow(x), same_key,
    isTRUE(all.equal(x, ref)), sum(x$outcome[kb] == ref$outcome[kb]), sum(kb), sum(x$confirmed[kp] == ref$confirmed[kp] & x$outcome[kp] == ref$outcome[kp]), sum(kp),
    sum(x$value[kb] == ref$value[kb]), sum(kb), sum(apply(x[!kb, ] == ref[!kb, ] | (is.na(x[!kb, ]) & is.na(ref[!kb, ])), 1, all), na.rm = TRUE), sum(!kb)) }
# B-2: were the modified RDS really used? compare b2/<v>/null summaries with S1.3 originals and with S1.3b res[[v]]
s13n <- list.files("/mnt/micron/geo_spatialtrans/S1.3/null", pattern = "[.]rds$")
chk <- mclapply(s13n, function(f) { o <- readRDS(file.path("/mnt/micron/geo_spatialtrans/S1.3/null", f)); b <- readRDS(file.path(B, "null", f))
  sapply(c("V0", "V1", "V2", "V3"), function(v) { x <- readRDS(file.path(B, "b2", v, "null", f))
    c(eq_s13 = isTRUE(all.equal(x$summary, o$summary, tolerance = 0)), eq_var = identical(unname(x$summary), unname(b$res[[v]]$summary)),
      cells_eq_var = if (is.null(x$cells)) NA else identical(x$cells$area, b$res[[v]]$cells$area)) }) }, mc.cores = 30)
A <- simplify2array(chk); L("B-2 inputs: summary identical to S1.3 original (per variant): %s", paste(apply(A["eq_s13", , ], 1, sum), collapse = "/"))
L("B-2 inputs: summary identical to S1.3b res[[v]]: %s; cells areas identical: %s", paste(apply(A["eq_var", , ], 1, sum), collapse = "/"), paste(apply(A["cells_eq_var", , ], 1, sum, na.rm = TRUE), collapse = "/"))
# B-2 margins: reproduce S1.3 check B (B1..B4) on each variant and compute minimum margin to thresholds
rs <- read.csv("results/R3/R3_roi_summary.csv"); prim <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
rs <- rs[rs$method == prim[rs$archetype], ]
key <- c("median_eq_r", "cv", "cv_loc", "median_ecc_T"); thr <- c(B1 = 0.05, B2 = 0.10, B3 = 0.10, B4 = 0.05)
mbs <- list()
for (v in c("V0", "V1", "V2", "V3")) {
  ns <- do.call(rbind, mclapply(s13n, function(f) { x <- readRDS(file.path(B, "b2", v, "null", f)); data.frame(archetype = x$archetype, roi_id = x$roi_id, model = x$model, t(x$summary[key])) }, mc.cores = 30))
  ag <- aggregate(ns[, key], ns[, c("archetype", "roi_id", "model")], median)
  mb <- merge(ag, rs[, c("archetype", "roi_id", key)], by = c("archetype", "roi_id"), suffixes = c("_sim", "_real"))
  mb$B1 <- mb$median_eq_r_sim / mb$median_eq_r_real - 1; mb$B2 <- mb$cv_sim / mb$cv_real - 1; mb$B3 <- mb$cv_loc_sim / mb$cv_loc_real - 1; mb$B4 <- mb$median_ecc_T_sim - mb$median_ecc_T_real
  mb <- mb[order(mb$archetype, mb$roi_id, mb$model), ]; mbs[[v]] <- mb }
for (v in c("V1", "V2", "V3")) { dmax <- sapply(names(thr), function(b) max(abs(mbs[[v]][[b]] - mbs$V0[[b]])))
  L("B-2 %s: max |ΔB| per metric %s", v, paste(names(dmax), signif(dmax, 3), collapse = ", ")) }
marg <- sapply(names(thr), function(b) min(abs(abs(mbs$V0[[b]]) - thr[[b]]))); L("B-2 V0 min margin | |B| - thr | over all ROI x model: %s", paste(names(marg), signif(marg, 3), collapse = ", "))
# B-2 order predictions & B-5 depend on cv/cv_loc (order) and cells (KS): margins
v0 <- mbs$V0; om <- c()
for (A in c("A3", "A5", "A6", "A4")) for (cc in c("cv", "cv_loc")) { g <- function(md) { d <- v0[v0$archetype == A & v0$model == md, ]; d[order(d$roi_id), paste0(cc, "_sim")] }
  om <- c(om, if (A == "A4") min(abs(c(g("RSArule") - g("RSA"), g("RSArule") - g("CSR")))) else min(abs(g("RSArule") - g("RSA")))) }
L("B-2 order predictions: min |G3 - G2 (or CSR)| margin %.3g (relative ~ %.3g)", min(om), min(om) / median(v0$cv_sim))
# ---- decision rule (recomputed)
el <- function(v) { e <- T[T$v == v, ]; gm <- do.call(rbind, lapply(G, function(o) o[[v]])); c2v <- C2[C2$variant == v, ]
  ut <- read.csv("results/S1.3b/S1.3b_test.csv"); u <- ut[ut$variant %in% c(v, "-") & ut$status != "INFO", ]
  c(C1 = sum(abs(e$ov) > e$tau_me | e$symd > e$tau_me) == 0, C2 = sum(e$pov) == 0, C3 = sum(e$cov_far) == 0 && sum(e$cov_anom) <= 50,
    C4 = sum(e$inv) == 0 && sum(e$emp) == 0 && all(e$n_out == e$n) && all(e$types_ok), C5 = sum(e$c5_over) == 0,
    C6 = sum(gm[, "cell_viol"]) + sum(gm[, "bin_viol"]) == 0, Cunit = all(u$status == "PASS"), CP2_posthoc = !any(c2v$fail_posthoc), CP2_silent = !any(c2v$silent),
    Cunit_n = nrow(u)) }
D <- sapply(c("V1", "V2", "V3"), el); print(D)
elig <- c("V1", "V3", "V2")[c("V1", "V3", "V2") %in% colnames(D)[apply(D[c("C1", "C2", "C3", "C4", "C5", "C6", "Cunit", "CP2_posthoc"), , drop = FALSE] == 1, 2, all)]]
L("eligible (order V1,V3,V2): %s; choice %s", paste(elig, collapse = ","), if (length(elig)) elig[1] else "P0")
elig2 <- c("V1", "V3", "V2")[c("V1", "V3", "V2") %in% colnames(D)[apply(D[c("C1", "C2", "C3", "C5", "C6", "Cunit", "CP2_posthoc"), , drop = FALSE] == 1, 2, all)]]
L("sensitivity: eligible ignoring C-4 invalid: %s", paste(elig2, collapse = ","))
