# results/S1.3/review/claims/01_recompute.R — revisione RA-claims: ricalcolo indipendente dei verdetti S1.3
# NON fa source di tools/S1.3_analysis.R. Legge solo uscite grezze e csv citati. Scrive solo in results/S1.3/review/claims/.
# Uso (dalla root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.3/review/claims/01_recompute.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(arrow))
RES <- "results/S1.3"; RAW <- "/mnt/micron/geo_spatialtrans/S1.3"; R3D <- "/mnt/micron/geo_spatialtrans/R3/real"
OD <- "results/S1.3/review/claims"; TOL <- 1e-9
prim <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
E2 <- function(v) formatC(v, format = "e", digits = 2)
PF <- function(x) ifelse(x, "PASS", "FAIL")
rmax <- function(a, b) if (!length(a)) 0 else max(abs(a - b) / pmax(abs(b), 1e-300))
R <- list(); add <- function(section, id, archetype, model, value, outcome, prediction = NA, confirmed = NA)
  R[[length(R) + 1]] <<- data.frame(section = section, id = id, archetype = archetype, model = model, value = as.character(value),
                                    outcome = outcome, prediction = as.character(prediction), confirmed = as.character(confirmed))
AUX <- list()   # numeri ausiliari per le verifiche del report

## ---- 1. test sintetici/avversari e mutanti -----------------------------------------------------------
tests <- list()
for (en in c("geos", "deldir")) {
  t <- read.csv(file.path(RES, sprintf("S1.3_test_%s_none.csv", en))); tests[[en]] <- t
  t <- t[t$status != "INFO", ]
  for (ck in unique(t$check)) { s <- t$status[t$check == ck]; add("C", ck, "-", en, sprintf("%d/%d", sum(s == "PASS"), length(s)), PF(all(s == "PASS"))) }
  AUX[[paste0("tot_assert_", en)]] <- sprintf("%d/%d", sum(t$status == "PASS"), nrow(t))
}
tg <- list(M1 = c("C-3", "C-5"), M2 = "C-1", M3 = "C-1", M4 = c("C-2", "C-7"), M5 = c("C-1", "C-4"), M6 = "C-3")   # dal testo pre-registrato
mut <- do.call(rbind, lapply(names(tg), function(m) { t <- read.csv(file.path(RES, sprintf("S1.3_test_geos_%s.csv", m)))
  data.frame(mutant = m, target_prereg = paste(tg[[m]], collapse = "/"), n_fail_target = sum(t$status == "FAIL" & t$check %in% tg[[m]]),
             n_fail_any = sum(t$status == "FAIL"), checks_failed = paste(unique(t$check[t$status == "FAIL"]), collapse = "/")) }))
write.csv(mut, file.path(OD, "rc_mutants.csv"), row.names = FALSE)
add("C", "mutanti", "-", "geos", sprintf("%d/6", sum(mut$n_fail_target > 0)), PF(all(mut$n_fail_target > 0)))

## ---- 2. ROI reali: C-5 ricalcolato dalle tabelle per cellula, C-1..3 ---------------------------------
rf <- sort(list.files(file.path(RAW, "real"), pattern = "\\.rds$", full.names = TRUE)); stopifnot(length(rf) == 40)
rr <- list(); c5x <- list(); realsum <- list(); realcells <- list(); arb <- list(); c10 <- list()
for (f in rf) {
  o <- readRDS(f)
  ref <- as.data.frame(read_parquet(file.path(R3D, sprintf("%s_%s_%s_cells.parquet", o$archetype, o$roi_id, o$method))))
  d <- o$geos$cells; stopifnot(nrow(d) == nrow(ref), all(d$idx == ref$idx))
  mono <- d$n_lost == 0 & d$n_gained == 0
  ii <- ref$interior & d$interior & mono
  i <- d$interior; a <- d$area[i]; al <- a * ref$lambda_loc[i]
  my_sum <- c(n = nrow(d), n_interior = sum(i), frac_interior = mean(i), mean_area = mean(a), median_area = median(a), cv = sd(a) / mean(a),
              median_eq_r = median(sqrt(a / pi)), median_ecc_T = median(d$ecc_T[i]), mean_nsides = mean(d$nsides[i]), var_nsides = var(d$nsides[i]),
              var_norm_area = var(a / mean(a)), cv_loc = sd(al) / mean(al))
  r3s <- read.csv("results/R3/R3_roi_summary.csv"); r3s <- r3s[r3s$archetype == o$archetype & r3s$roi_id == o$roi_id & r3s$method == o$method, ]
  c5x[[length(c5x) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, method = o$method,
    my_interior_identical = identical(d$interior, ref$interior), my_n_interior_diff = sum(d$interior != ref$interior),
    my_rel_area_interior = rmax(d$area[ref$interior], ref$area[ref$interior]), my_rel_area_mono = rmax(d$area[mono], ref$area[mono]),
    my_rel_area_interior_mono = rmax(d$area[ref$interior & mono], ref$area[ref$interior & mono]),
    n_frag_cells = sum(!mono), n_interior_frag = sum(ref$interior & !mono), n_interior_diff_mono = sum((d$interior != ref$interior) & mono),
    my_rel_sum_area = abs(sum(d$area) - sum(ref$area)) / sum(ref$area),
    nsides_eq_int_mono = mean(d$nsides[ii] == ref$nsides[ii]), max_abs_decc_int_mono = max(abs(d$ecc_T[ii] - ref$ecc_T[ii])),
    my_rel_summary_vs_R3 = rmax(my_sum, unlist(r3s[names(my_sum)])), my_summary_vs_stored = rmax(my_sum, o$geos$summary[names(my_sum)]),
    stored_rel_summary = o$geos$c5[["rel_summary"]])
  for (en in c("geos", "deldir")) { e <- o[[en]]
    rr[[length(rr) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, method = o$method, engine = en, n = o$n, t(e$c123), t(e$c5)) }
  if (o$method == prim[[o$archetype]]) { realsum[[length(realsum) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, t(o$geos$summary))
    realcells[[length(realcells) + 1]] <- data.frame(archetype = o$archetype, a = a / mean(a), nsides = d$nsides[i]) }
  if (nrow(o$arb)) arb[[length(arb) + 1]] <- data.frame(set = "real", archetype = o$archetype, roi_id = o$roi_id, model = o$method, o$arb)
  c10[[length(c10) + 1]] <- data.frame(set = "real", archetype = o$archetype, roi_id = o$roi_id, model = o$method, t(o$c10a[1:4]))
}
rr <- do.call(rbind, rr); c5x <- do.call(rbind, c5x); realsum <- do.call(rbind, realsum); realcells <- do.call(rbind, realcells)
rr$c5_text <- rr$interior_identical == 1 & rr$rel_area_interior <= TOL & rr$rel_area_mono <= TOL & rr$rel_sum_area <= TOL   # testo pre-registrato
rr$c5_coded <- rr$c5_text & rr$rel_summary <= TOL                                                                             # + identita' del riassunto
rr$c123 <- rr$c1 <= TOL & rr$c2_overlap <= TOL & rr$c2_symdiff <= TOL & rr$c3_valid == 1 & rr$c3_gen_in_own == 1 & rr$n_repaired == 0
rr$c123_strict <- rr$c123 & rr$n_multipart == 0
write.csv(rr, file.path(OD, "rc_roi_real.csv"), row.names = FALSE); write.csv(c5x, file.path(OD, "rc_c5_from_cells.csv"), row.names = FALSE)
for (en in c("geos", "deldir")) { k <- rr$engine == en
  add("C", "C-5", "A1-A6", en, sprintf("%d/%d; max rel area interne %s", sum(rr$c5_coded[k]), sum(k), E2(max(rr$rel_area_interior[k]))), PF(all(rr$c5_coded[k])))
  add("C", "C-1..3 ROI", "A1-A6", en, sprintf("%d/%d; max C-1 %s", sum(rr$c123[k]), sum(k), E2(max(rr$c1[k]))), PF(all(rr$c123[k]))) }

## ---- 3. nulli -------------------------------------------------------------------------------------------
nf <- sort(list.files(file.path(RAW, "null"), pattern = "\\.rds$", full.names = TRUE)); stopifnot(length(nf) == 1800)
ns <- list(); dn <- list(); simcells <- list(); recheck <- list()
for (f in nf) { o <- readRDS(f)
  ns[[length(ns) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep, n_failed = o$n_failed, t(o$summary), t(o$c123))
  if (!is.null(o$c123_deldir)) dn[[length(dn) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep, t(o$c123_deldir))
  if (!is.null(o$c10a)) c10[[length(c10) + 1]] <- data.frame(set = "null", archetype = o$archetype, roi_id = o$roi_id, model = o$model, t(o$c10a[1:4]))
  if (!is.null(o$arb) && nrow(o$arb)) arb[[length(arb) + 1]] <- data.frame(set = "null", archetype = o$archetype, roi_id = o$roi_id, model = o$model, o$arb)
  if (!is.null(o$cells)) { d <- o$cells; i <- d$interior; a <- d$area[i]
    simcells[[length(simcells) + 1]] <- data.frame(archetype = o$archetype, model = o$model, rep = o$rep, a = a / mean(a), nsides = d$nsides[i])
    ms <- c(n = nrow(d), n_interior = sum(i), mean_area = mean(a), median_area = median(a), cv = sd(a) / mean(a), median_eq_r = median(sqrt(a / pi)),
            median_ecc_T = median(d$ecc_T[i]), mean_nsides = mean(d$nsides[i]), var_nsides = var(d$nsides[i]), var_norm_area = var(a / mean(a)))
    recheck[[length(recheck) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep,
      rel_vs_stored = rmax(ms, o$summary[names(ms)]), nsides_hist_ok = identical(as.integer(tabulate(pmin(d$nsides[i], 15), 15)), as.integer(o$nsides))) }
}
ns <- do.call(rbind, ns); dn <- do.call(rbind, dn); simcells <- do.call(rbind, simcells); recheck <- do.call(rbind, recheck)
write.csv(recheck, file.path(OD, "rc_null_summary_from_cells.csv"), row.names = FALSE)
cols <- c("n", "n_interior", "frac_interior", "mean_area", "median_area", "cv", "median_eq_r", "median_ecc_T", "mean_nsides", "var_nsides", "var_norm_area", "cv_loc")
r3n <- read.csv("results/R3/R3_null_summary.csv")
m6 <- merge(ns[ns$model %in% c("CSR", "RSA"), c("archetype", "roi_id", "model", "rep", "n_failed", "n_fragments", cols)], r3n[, c("archetype", "roi_id", "model", "rep", "n_failed", cols)],
            by = c("archetype", "roi_id", "model", "rep"), suffixes = c("", ".r3"))
rel <- sapply(cols, function(cc) abs(m6[[cc]] - m6[[paste0(cc, ".r3")]]) / pmax(abs(m6[[paste0(cc, ".r3")]]), 1e-300))
m6$rel <- apply(rel, 1, max); m6$worst <- cols[apply(rel, 1, which.max)]; m6$ident <- m6$rel <= TOL & m6$n_failed == m6$n_failed.r3
m6$rel_no_interior_metrics <- apply(rel[, c("n"), drop = FALSE], 1, max)
write.csv(m6[, c("archetype", "roi_id", "model", "rep", "n_fragments", "rel", "worst", "ident")], file.path(OD, "rc_c6.csv"), row.names = FALSE)
add("C", "C-6", "A1-A6", "geos", sprintf("%d/%d (attese 1200); max rel %s", sum(m6$ident), nrow(m6), E2(max(m6$rel))), PF(nrow(m6) == 1200 && all(m6$rel <= TOL)))
AUX$c6_nofrag <- sprintf("%d/%d", sum(m6$ident[m6$n_fragments == 0]), sum(m6$n_fragments == 0))
AUX$c6_frag <- sprintf("%d/%d", sum(m6$ident[m6$n_fragments > 0]), sum(m6$n_fragments > 0))
AUX$c6_nofrag_maxrel_nonident <- E2(max(c(0, m6$rel[m6$n_fragments == 0 & !m6$ident])))
ns$c123 <- ns$c1 <= TOL & ns$c2_overlap <= TOL & ns$c2_symdiff <= TOL & ns$c3_valid == 1 & ns$c3_gen_in_own == 1 & ns$n_repaired == 0
dn$c123 <- dn$c1 <= TOL & dn$c2_overlap <= TOL & dn$c2_symdiff <= TOL & dn$c3_valid == 1 & dn$c3_gen_in_own == 1 & dn$n_repaired == 0
add("C", "C-1..3 nulli", "A1-A6", "geos", sprintf("%d/%d", sum(ns$c123), nrow(ns)), PF(nrow(ns) == 1800 && all(ns$c123)))
add("C", "C-1..3 nulli", "A1-A6", "deldir", sprintf("%d/%d", sum(dn$c123), nrow(dn)), PF(nrow(dn) == 90 && all(dn$c123)))
AUX$null_geos_fail_criteria <- paste(sprintf("%s=%d", c("c1", "c2_overlap", "c2_symdiff", "c3_valid", "gen_in_own", "repaired", "multipart>0"),
  c(sum(ns$c1 > TOL), sum(ns$c2_overlap > TOL), sum(ns$c2_symdiff > TOL), sum(ns$c3_valid != 1), sum(ns$c3_gen_in_own != 1), sum(ns$n_repaired > 0), sum(ns$n_multipart > 0))), collapse = " ")
AUX$null_deldir_fail_criteria <- paste(sprintf("%s=%d", c("c1", "c2_overlap", "c2_symdiff", "c3_valid", "gen_in_own", "repaired", "multipart>0"),
  c(sum(dn$c1 > TOL), sum(dn$c2_overlap > TOL), sum(dn$c2_symdiff > TOL), sum(dn$c3_valid != 1), sum(dn$c3_gen_in_own != 1), sum(dn$n_repaired > 0), sum(dn$n_multipart > 0))), collapse = " ")
AUX$null_geos_strict_with_multipart <- sprintf("%d/%d", sum(ns$c123 & ns$n_multipart == 0), nrow(ns))
AUX$null_deldir_strict_with_multipart <- sprintf("%d/%d", sum(dn$c123 & dn$n_multipart == 0), nrow(dn))
AUX$null_n_failed_gt0 <- sum(ns$n_failed > 0)
AUX$null_reps_per_roi <- paste(range(table(paste(ns$archetype, ns$roi_id, ns$model))), collapse = "-")
# gli stessi 90 casi replica 1 per geos, per un confronto di tasso alla pari
g1 <- merge(ns[ns$rep == 1, c("archetype", "roi_id", "model", "c123")], dn[, c("archetype", "roi_id", "model", "c123")], by = c("archetype", "roi_id", "model"), suffixes = c("_geos", "_deldir"))
AUX$rep1_c123_geos_vs_deldir <- sprintf("geos %d/%d, deldir %d/%d", sum(g1$c123_geos), nrow(g1), sum(g1$c123_deldir), nrow(g1))

## ---- C-10a e arbitro ------------------------------------------------------------------------------------
syn <- read.csv(file.path(RES, "S1.3_c10a_syn.csv"))
c10 <- rbind(do.call(rbind, c10), data.frame(set = syn$set, archetype = syn$archetype, roi_id = syn$roi_id, model = syn$model, rel_area = syn$rel_area,
             interior_identical = as.numeric(as.logical(syn$interior_identical)), frag_identical = as.numeric(as.logical(syn$frag_identical)),
             nsides_interior_identical = as.numeric(as.logical(syn$nsides_interior_identical))))
c10$ok <- c10$rel_area <= TOL & c10$interior_identical == 1 & c10$frag_identical == 1 & c10$nsides_interior_identical == 1
arb <- do.call(rbind, arb)
arb$my_geos_ok <- abs(arb$a_geos - arb$a_hp) / arb$a_hp <= 1e-9 & arb$ns_geos == arb$ns_hp
arb$my_deldir_ok <- abs(arb$a_deldir - arb$a_hp) / arb$a_hp <= 1e-9 & arb$ns_deldir == arb$ns_hp
write.csv(arb, file.path(OD, "rc_arbitration.csv"), row.names = FALSE)
nd <- nrow(arb) + sum(syn$n_discord); okg <- sum(arb$geos_ok) + sum(syn$n_geos_ok); okd <- sum(arb$deldir_ok) + sum(syn$n_deldir_ok)
add("C10", "arbitro", "A1-A6", "semipiani", sprintf("%d: %d / %d (limitati %d)", nd, okg, okd, sum(arb$bounded)), "INFO")
add("C10", "C-10a", "A1-A6", "geos vs deldir", sprintf("%d/%d; max rel area %s", sum(c10$ok), nrow(c10), E2(max(c10$rel_area))), PF(all(c10$ok)), "PASS", all(c10$ok))
AUX$c10a_sets <- sprintf("real %d, null %d, syn %d", sum(c10$set == "real"), sum(c10$set == "null"), sum(c10$set == "syn"))
AUX$c10a_fail_by_set <- paste(tapply(!c10$ok, c10$set, sum), collapse = "/")
AUX$arb_min_edge <- paste(E2(range(arb$min_edge)), collapse = " .. ")
AUX$arb_myflags_agree <- sprintf("geos %d/%d, deldir %d/%d", sum(arb$my_geos_ok == arb$geos_ok), nrow(arb), sum(arb$my_deldir_ok == arb$deldir_ok), nrow(arb))
AUX$arb_geos_wrong <- paste(apply(arb[!arb$geos_ok, c("set", "archetype", "roi_id", "model", "cell", "a_geos", "a_hp", "ns_geos", "ns_hp", "ns_deldir")], 1, paste, collapse = ":"), collapse = " | ")
AUX$arb_area_only_wrong_deldir <- sum(!arb$deldir_ok & abs(arb$a_deldir - arb$a_hp) / arb$a_hp > 1e-9)

## ---- 4. check B -----------------------------------------------------------------------------------------
rs <- read.csv("results/R3/R3_roi_summary.csv"); rs <- rs[rs$method == prim[rs$archetype], ]; stopifnot(nrow(rs) == 30)
key <- c("median_eq_r", "cv", "cv_loc", "median_ecc_T")
agg <- aggregate(ns[, key], ns[, c("archetype", "roi_id", "model")], median)
mkB <- function(ref) { m <- merge(agg, ref[, c("archetype", "roi_id", key)], by = c("archetype", "roi_id"), suffixes = c("_sim", "_real"))
  m$B1 <- m$median_eq_r_sim / m$median_eq_r_real - 1; m$B2 <- m$cv_sim / m$cv_real - 1; m$B3 <- m$cv_loc_sim / m$cv_loc_real - 1
  m$B4 <- m$median_ecc_T_sim - m$median_ecc_T_real
  m$p1 <- abs(m$B1) <= 0.05; m$p2 <- abs(m$B2) <= 0.10; m$p3 <- abs(m$B3) <= 0.10; m$p4 <- abs(m$B4) < 0.05; m }
mb <- mkB(rs); ms <- mkB(realsum)
write.csv(mb, file.path(OD, "rc_checkB_roi.csv"), row.names = FALSE)
# previsioni trascritte dal TESTO della pre-registrazione (non dal codice)
P <- list(RSA = list(B1 = c("FAIL", "FAIL", "FAIL", "FAIL", "FAIL", "PASS"), B2 = c(rep("FAIL", 5), "PASS"), B3 = c(rep("FAIL", 5), "PASS"), B4 = rep("PASS", 6)),
          RSArule = list(B1 = c("FAIL", "FAIL", NA, NA, NA, NA), B2 = c(NA, NA, NA, NA, NA, "FAIL"), B3 = c(NA, NA, NA, NA, NA, "FAIL"), B4 = c(NA, NA, NA, NA, "FAIL", "FAIL")))
pr3 <- read.csv(file.path(RES, "S1.3_predictions_from_R3.csv")); pr3 <- pr3[pr3$model == "CSR", ]
G <- c(CSR = "G1", RSA = "G2", RSArule = "G3"); g1pred <- list()
for (md in c("CSR", "RSA", "RSArule")) for (A in paste0("A", 1:6)) for (b in c("B1", "B2", "B3", "B4")) {
  d <- mb[mb$archetype == A & mb$model == md, ]; stopifnot(nrow(d) == 5); np <- sum(d[[sub("B", "p", b)]]); out <- PF(np >= 4)
  p <- if (md == "CSR") NA else P[[md]][[b]][as.integer(substr(A, 2, 2))]
  add("B", b, A, G[[md]], sprintf("%d/5 ROI; %+.3f..%+.3f", np, min(d[[b]]), max(d[[b]])), out, p, if (is.na(p)) NA else p == out)
  if (md == "CSR") { pn <- pr3[pr3$archetype == A, paste0(b, "_n_pass")]
    g1pred[[length(g1pred) + 1]] <- data.frame(archetype = A, metric = b, pred_n_pass_R3 = pn, pred = PF(pn >= 4), obs_n_pass = np, obs = out, confirmed = PF(pn >= 4) == out) }
}
g1pred <- do.call(rbind, g1pred); write.csv(g1pred, file.path(OD, "rc_G1_predictions_uncoded.csv"), row.names = FALSE)
for (md in c("RSA", "RSArule")) for (A in paste0("A", 1:6)) for (b in c("B1", "B2", "B3", "B4")) {
  d <- ms[ms$archetype == A & ms$model == md, ]; d0 <- mb[mb$archetype == A & mb$model == md, ]
  o1 <- PF(sum(d[[sub("B", "p", b)]]) >= 4); o0 <- PF(sum(d0[[sub("B", "p", b)]]) >= 4)
  add("Bsens", b, A, G[[md]], sprintf("%d/5 ROI", sum(d[[sub("B", "p", b)]])), o1, o0, o1 == o0) }
for (A in c("A3", "A5", "A6", "A4")) for (cc in c("cv", "cv_loc")) {
  v <- function(md) { d <- mb[mb$archetype == A & mb$model == md, ]; d[order(d$roi_id), paste0(cc, "_sim")] }
  v1 <- v("CSR"); v2 <- v("RSA"); v3 <- v("RSArule")
  ok <- if (A == "A4") sum(v3 > pmin(v1, v2) & v3 < pmax(v1, v2)) else sum(v3 < v2)
  add("B", paste0("ordine ", cc), A, "G3", sprintf("%d/5 ROI", ok), PF(ok >= 4), "PASS", ok >= 4) }
ksd <- function(a, b) suppressWarnings(unname(ks.test(a, b)$statistic))
ksn <- list()
for (A in paste0("A", 1:6)) for (md in c("CSR", "RSA", "RSArule")) {
  re <- realcells[realcells$archetype == A, ]; s1 <- simcells[simcells$archetype == A & simcells$model == md & simcells$rep == 1, ]
  s2 <- simcells[simcells$archetype == A & simcells$model == md & simcells$rep == 2, ]
  ksn[[length(ksn) + 1]] <- data.frame(archetype = A, model = md, D = ksd(re$a, s1$a), D_noise = ksd(s1$a, s2$a), var_nsides_real = var(re$nsides), var_nsides_sim = var(s1$nsides))
  add("B", "B-5", A, G[[md]], sprintf("%.3f (%.3f)", ksd(re$a, s1$a), ksd(s1$a, s2$a)), "INFO") }
write.csv(do.call(rbind, ksn), file.path(OD, "rc_ks_nsides.csv"), row.names = FALSE)

## ---- 5. controprove e descrittivi ----------------------------------------------------------------------
cp1 <- read.csv(file.path(RES, "S1.3_cp1_regions.csv")); k <- cp1$n_cells >= 3 & cp1$n_cells <= 1000
fd <- mean(cp1$uncovered_design[k] > 0.01); npk <- sum(abs(cp1$uncovered_pkg[cp1$n_cells > 0]) > TOL)
add("CP", "CP-1", "sintetici", "testo del design", sprintf("%.3f (n=%d); mediana %.4f", fd, sum(k), median(cp1$uncovered_design[k])), PF(fd >= 0.5), "PASS", fd >= 0.5)
add("CP", "CP-1", "sintetici", "D-S1.3.1", npk, PF(npk == 0), "PASS", npk == 0)
AUX$cp1_seeds <- length(unique(cp1$seed)); AUX$cp1_inputs <- paste(unique(cp1$input), collapse = ",")
AUX$cp1_ncells_quantiles <- paste(round(quantile(cp1$n_cells, c(0, .25, .5, .75, .9, 1))), collapse = "/")
AUX$cp1_regions_total <- nrow(cp1); AUX$cp1_regions_gt1000 <- sum(cp1$n_cells > 1000); AUX$cp1_regions_lt3 <- sum(cp1$n_cells < 3)
AUX$cp1_median_all_ge3 <- sprintf("%.4f", median(cp1$uncovered_design[cp1$n_cells >= 3]))
AUX$cp1_area_quantiles_k <- paste(round(quantile(cp1$area_um2[k], c(0, .5, 1))), collapse = "/")
cpk <- function(fn) { x <- read.csv(file.path(RES, fn)); a <- aggregate(median_ecc_T ~ k, x, median); a$delta <- a$median_ecc_T - a$median_ecc_T[a$k == 1]
  list(a = a, mono = all(diff(a$median_ecc_T) > 0), kmin = suppressWarnings(min(a$k[a$delta >= 0.05])), nrep = paste(range(table(x$k)), collapse = "-")) }
c2a <- cpk("S1.3_cp2_anisotropy.csv"); c2b <- cpk("S1.3_cp2b_posthoc_regular.csv")
add("CP", "CP-2", "A6 (densita')", "RSA affine", paste(sprintf("%.3f", c2a$a$median_ecc_T), collapse = " / "), PF(c2a$mono), "PASS", c2a$mono)
add("CP", "CP-2", "A6 (densita')", "RSA affine", c2a$kmin, PF(is.finite(c2a$kmin) && c2a$kmin <= 1.5), "PASS", is.finite(c2a$kmin) && c2a$kmin <= 1.5)
add("CP", "CP-2b (post hoc)", "A6 (densita')", "RSA regolare affine", paste(sprintf("%.3f", c2b$a$median_ecc_T), collapse = " / "), PF(c2b$mono), "PASS", c2b$mono)
add("CP", "CP-2b (post hoc)", "A6 (densita')", "RSA regolare affine", c2b$kmin, PF(is.finite(c2b$kmin) && c2b$kmin <= 1.5), "PASS", is.finite(c2b$kmin) && c2b$kmin <= 1.5)
AUX$cp2_delta <- paste(sprintf("%.4f", c2a$a$delta), collapse = "/"); AUX$cp2b_delta <- paste(sprintf("%.4f", c2b$a$delta), collapse = "/")
AUX$cp2_reps <- c2a$nrep; AUX$cp2b_reps <- c2b$nrep
cf <- sort(list.files(file.path(RAW, "cp3"), pattern = "\\.rds$", full.names = TRUE))
c3 <- do.call(rbind, lapply(cf, function(f) { o <- readRDS(f); data.frame(archetype = o$archetype, roi_id = o$roi_id, ovc = o$ov_contained$n, ovm = o$ov_m4$n,
  bad = o$ov_m4$n_bad, nnc = o$n_clipped_nonconvex, g13 = o$gap_frac[[1]], g23 = o$gap_frac[[2]], g1 = o$gap_frac[[3]]) }))
write.csv(c3, file.path(OD, "rc_cp3.csv"), row.names = FALSE)
add("CP", "CP-3", "A1-A6", "Chaikin senza contenimento", sprintf("%d (%d/%d)", sum(c3$ovm), sum(c3$ovm > 0), nrow(c3)), PF(sum(c3$ovm) > 0), "PASS", sum(c3$ovm) > 0)
add("CP", "CP-3", "A1-A6", "Chaikin senza contenimento", sum(c3$bad), PF(sum(c3$bad) == 0), if (sum(c3$ovm) == 0) "PASS (vacua: 0 coppie)" else "PASS", sum(c3$bad) == 0)
add("CP", "CP-3", "A1-A6", "D-S1.3.3", sum(c3$ovc), PF(sum(c3$ovc) == 0), "PASS", sum(c3$ovc) == 0)
gm <- aggregate(g1 ~ archetype, c3, median); AUX$cp3_gap_cs1_range <- sprintf("%.4f..%.4f", min(gm$g1), max(gm$g1))
AUX$cp3_gap_monotone <- all(c3$g13 < c3$g23 & c3$g23 < c3$g1)
d1 <- read.csv(file.path(RES, "S1.3_d1_boundary.csv"))
add("D", "D-1", "A1 (densita')", "RSA regola", sprintf("%.3f", median(d1$ratio_band)), "INFO", ">= 1.1", median(d1$ratio_band) >= 1.1)
AUX$d1 <- sprintf("n=%d median %.3f min %.3f max %.3f; pooled ratio mean(band_r2)/mean(band_r1)=%.3f; seeds>=1.1: %d", nrow(d1), median(d1$ratio_band), min(d1$ratio_band), max(d1$ratio_band),
                  mean(d1$band_r2) / mean(d1$band_r1), sum(d1$ratio_band >= 1.1))
w <- cp1$n_cells > 0
add("D", "D-3", "sintetici", "-", sprintf("%.2f", suppressWarnings(cor(cp1$clipped_frac[w], log(cp1$area_um2[w]), method = "spearman"))), "INFO")
pf <- read.csv(file.path(RES, "S1.3_perf.csv"))
sl <- sapply(c("geos", "deldir"), function(b) { d <- pf[pf$backend == b & pf$n <= 160000, ]; unname(coef(lm(log(elapsed_s) ~ log(n), d))[2]) })
for (b in c("geos", "deldir")) { d <- pf[pf$backend == b & pf$n == 449536, ]
  ok <- d$status == "ok" && d$elapsed_s <= 150 && d$max_rss_mb <= 4096
  add("C", "C-perf", "-", b, sprintf("%s s, %s MB (%s)", d$elapsed_s, d$max_rss_mb, d$status), PF(ok))
  add("C10", "C-10c", "-", b, sprintf("%.2f", sl[[b]]), "INFO", if (b == "geos") "≈ 1" else "≈ 2") }
rat <- merge(pf[pf$backend == "deldir", c("n", "elapsed_s", "max_rss_mb")], pf[pf$backend == "geos", c("n", "elapsed_s", "max_rss_mb")], by = "n", suffixes = c("_deldir", "_geos"))
rat$speedup <- rat$elapsed_s_deldir / rat$elapsed_s_geos; rat$mem_ratio <- rat$max_rss_mb_deldir / rat$max_rss_mb_geos
write.csv(rat, file.path(OD, "rc_perf_ratio.csv"), row.names = FALSE)
AUX$speedup_by_n <- paste(sprintf("%d:%.1fx", rat$n, rat$speedup), collapse = " ")

## ---- scelta del motore (regola dell'addendum) ---------------------------------------------------------
rob <- function(en) { t <- tests[[en]][tests[[en]]$check == "ROB", ]
  st <- if (nrow(t)) sum(as.integer(unlist(strsplit(t$value, "/")))) else 0
  sr <- with(rr[rr$engine == en, ], sum(n_snapped + n_multipart + n_repaired))
  sn <- if (en == "geos") with(ns, sum(n_snapped + n_multipart + n_repaired)) else with(dn, sum(n_snapped + n_multipart + n_repaired))
  c(test = st, real = sr, null = sn) }
K <- sapply(c("geos", "deldir"), function(en) { t <- tests[[en]]; t <- t[t$status != "INFO", ]
  c(K1 = all(t$status == "PASS"), K2 = all(rr$c5_coded[rr$engine == en] & rr$c123[rr$engine == en]),
    K2_text = all(rr$c5_text[rr$engine == en] & rr$c123[rr$engine == en]),
    K3 = if (en == "geos") all(ns$c123) else all(dn$c123), Karb = (if (en == "geos") okg else okd) == nd) })
elig <- colnames(K)[K["K1", ] & K["K2", ] & K["K3", ]]
add("C10", "scelta", "-", paste(elig, collapse = "+"), if (length(elig)) elig[1] else NA, if (length(elig)) "INFO" else "FAIL", "geos", if (length(elig) == 1) elig == "geos" else FALSE)
eng <- data.frame(engine = colnames(K), t(K), rob_test = sapply(colnames(K), function(e) rob(e)[["test"]]), rob_real = sapply(colnames(K), function(e) rob(e)[["real"]]),
                  rob_null = sapply(colnames(K), function(e) rob(e)[["null"]]), n_null_tess = c(nrow(ns), nrow(dn)),
                  snapped_real = sapply(colnames(K), function(e) sum(rr$n_snapped[rr$engine == e])), multipart_real = sapply(colnames(K), function(e) sum(rr$n_multipart[rr$engine == e])),
                  elapsed_roi_range = sapply(colnames(K), function(e) paste(range(rr$elapsed_s[rr$engine == e]), collapse = "-")))
write.csv(eng, file.path(OD, "rc_engine_choice.csv"), row.names = FALSE)
# scomposizione di K-2 per criterio
k2 <- do.call(rbind, lapply(c("geos", "deldir"), function(en) { d <- rr[rr$engine == en, ]
  data.frame(engine = en, fail_interior_flags = sum(d$interior_identical != 1), fail_rel_area_interior = sum(d$rel_area_interior > TOL),
             fail_rel_area_mono = sum(d$rel_area_mono > TOL), max_rel_area_mono = max(d$rel_area_mono), fail_rel_sum = sum(d$rel_sum_area > TOL),
             fail_rel_summary = sum(d$rel_summary > TOL), fail_c1 = sum(d$c1 > TOL), fail_c2_overlap = sum(d$c2_overlap > TOL), fail_c2_symdiff = sum(d$c2_symdiff > TOL),
             fail_c3 = sum(d$c3_valid != 1 | d$c3_gen_in_own != 1), multipart_roi = sum(d$n_multipart > 0),
             c5_pass_text = sum(d$c5_text), c5_pass_coded = sum(d$c5_coded), c5_pass_if_mono_only = sum(d$rel_area_mono <= TOL & d$rel_sum_area <= TOL),
             roi_c123_fail = paste(paste(d$archetype, d$roi_id, d$method)[!d$c123], collapse = "; ")) }))
write.csv(k2, file.path(OD, "rc_K2_decomposition.csv"), row.names = FALSE)
cv <- read.csv(file.path(RES, "S1.3_c2_coverage.csv"))
for (e in c("geos", "deldir")) { d <- cv[cv$engine == e, ]
  add("C", "C-2 copertura (post hoc)", "A1-A6 + I4", e, sprintf("%d/%d/%d; %d coppie; %.2e um2 (%d casi)", sum(d$cov0), sum(d$cov1), sum(d$cov2), sum(d$n_pairs), sum(d$ov_area), nrow(d)),
      PF(sum(d$cov0) == 0 && sum(d$cov2) == 0)) }
AUX$coverage_sets <- paste(capture.output(print(table(cv$set, cv$engine))), collapse = " ; ")
AUX$coverage_pts_per_case <- paste(unique(cv$n_pts), collapse = ",")
pc <- read.csv(file.path(RES, "S1.3_c2_overlap_pairs_check.csv")); AUX$pairs_check_rows <- nrow(pc); AUX$pairs_check_le9e13 <- sum(pc$area <= 9e-13)
AUX$geos_cov_pairs_by_set <- paste(tapply(cv$n_pairs[cv$engine == "geos"], cv$set[cv$engine == "geos"], sum), collapse = "/")
a5 <- cv[cv$archetype == "A5" & cv$roi_id == "r5", ]; AUX$a5r5_cov <- paste(apply(a5[, c("engine", "case", "c2_union")], 1, paste, collapse = ":"), collapse = " | ")

MY <- do.call(rbind, R); write.csv(MY, file.path(OD, "rc_my_verdicts.csv"), row.names = FALSE)
writeLines(paste(names(AUX), unlist(lapply(AUX, function(x) paste(x, collapse = " "))), sep = " = "), file.path(OD, "rc_aux.txt"))
## ---- confronto con S1.3_verdicts.csv ---------------------------------------------------------------------
V <- read.csv(file.path(RES, "S1.3_verdicts.csv"), check.names = FALSE, colClasses = "character")
kk <- function(d) { k <- paste(d$section, d$id, d$archetype, d$model, sep = "|"); paste(k, ave(seq_along(k), k, FUN = seq_along), sep = "#") }
V$key <- kk(V); MY$key <- kk(MY)
cmp <- merge(V[, c("key", "value", "outcome", "prediction", "confirmed")], MY[, c("key", "value", "outcome", "prediction", "confirmed")], by = "key", all = TRUE, suffixes = c("_orig", "_mine"))
nz <- function(x) ifelse(is.na(x) | x == "NA", "", x)
cmp$same_value <- nz(cmp$value_orig) == nz(cmp$value_mine); cmp$same_outcome <- nz(cmp$outcome_orig) == nz(cmp$outcome_mine)
cmp$same_pred <- nz(cmp$prediction_orig) == nz(cmp$prediction_mine); cmp$same_conf <- nz(cmp$confirmed_orig) == nz(cmp$confirmed_mine)
cmp$identical <- cmp$same_value & cmp$same_outcome & cmp$same_pred & cmp$same_conf
write.csv(cmp, file.path(OD, "rc_compare_verdicts.csv"), row.names = FALSE)
cat("rows orig", nrow(V), "mine", nrow(MY), "merged", nrow(cmp), "\n")
cat("identical all fields:", sum(cmp$identical, na.rm = TRUE), " same outcome:", sum(cmp$same_outcome, na.rm = TRUE), "\n")
print(cmp[!cmp$identical | is.na(cmp$identical), c("key", "value_orig", "value_mine", "outcome_orig", "outcome_mine", "prediction_orig", "prediction_mine", "confirmed_orig", "confirmed_mine")], row.names = FALSE)
cat(readLines(file.path(OD, "rc_aux.txt")), sep = "\n")
print(eng, row.names = FALSE); print(k2[, -ncol(k2)], row.names = FALSE); print(k2$roi_c123_fail)
print(g1pred, row.names = FALSE)
print(mut, row.names = FALSE)
