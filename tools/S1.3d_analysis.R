# tools/S1.3d_analysis.R — S1.3d (BL-084): aggregazione, verdetti e regola di decisione (pre-registrata, ee6208b).
# Uso (root del repo): Rscript --vanilla tools/S1.3d_analysis.R [outdir] [--holdout]
# Legge: res/*.rds (candidati vs verità), coverage/*.json, geom/*_meta.json, checks/arbiter_check.json, res/mutants.rds, logs/perf_*.log.
# Scrive in results/S1.3d/: S1.3d_o1_long.csv, S1.3d_o2_long.csv, S1.3d_perf.csv, S1.3d_eligibility.csv, S1.3d_decision.csv, S1.3d_verdicts.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(jsonlite); library(arrow) })
a <- commandArgs(TRUE); OUTD <- if (length(a) && !startsWith(a[1], "--")) a[1] else "/mnt/micron/geo_spatialtrans/S1.3d"; HOLD <- "--holdout" %in% a
RES <- file.path(OUTD, "res"); GEOM <- file.path(OUTD, "geom"); COV <- file.path(OUTD, "coverage"); OUTR <- Sys.getenv("S13D_OUTR", "results/S1.3d"); dir.create(OUTR, showWarnings = FALSE)
G_GRID <- 1e-6; GUAR <- list(O1 = c("F1", "F5", "F4"), O2 = c("F1", "F5"), O3 = c("F1", "F5"))
CANDS <- list(O1 = c("V0", "F1", "F2", "F3", "F4", "F5"), O2 = c("V0", "F1", "F2", "F3", "F4", "F5", "F6"), O3 = c("V0", "F1", "F2", "F3", "F4", "F5", "F6"))
set_of <- function(k) suppressWarnings(ifelse(grepl("^U1_", k), "U1", ifelse(grepl("^U2_", k), "U2", ifelse(grepl("^U3_", k), "U3", ifelse(grepl("^U4_", k), "U4", ifelse(grepl("^U5_", k), "U5",
          ifelse(grepl("^CP1T_", k), "CP1T", ifelse(grepl("^CP1P_", k), "CP1P", ifelse(grepl("^S_", k), "S", ifelse(as.integer(sub(".*_", "", k)) >= 61, "H", "N"))))))))))
rbind_fill <- function(L) { L <- Filter(function(d) !is.null(d) && nrow(d) > 0, L); nm <- unique(unlist(lapply(L, names)))
  do.call(rbind, lapply(L, function(d) { for (k in setdiff(nm, names(d))) d[[k]] <- NA; d[, nm, drop = FALSE] })) }

# ---- lettura -------------------------------------------------------------------------------------------------------
fs <- list.files(RES, pattern = "_o12?\\.rds$", full.names = TRUE)
L <- lapply(fs, readRDS)
O1 <- rbind_fill(lapply(L, function(o) { d <- o$o1; if (is.null(d)) return(NULL); d$key <- o$key; d }))
O2 <- rbind_fill(lapply(L, function(o) { d <- o$o2; if (!length(d)) return(NULL); d$key <- o$key; d }))
C5 <- rbind_fill(lapply(L, function(o) if (is.null(o$c5) || !is.null(o$c5$error)) { if (!is.null(o$c5)) data.frame(key = o$key, error = o$c5$error) } else
  data.frame(key = o$key, n_over = o$c5$n_over, max_ratio = o$c5$max_ratio, same_lost = o$c5$same_lost, same_gained = o$c5$same_gained,
             n_tie = o$c5$info$n_tie, min_margin = o$c5$info$min_margin, n_ambiguous = o$c5$info$n_ambiguous, n_invalid = o$c5$info$n_invalid, t_build = o$c5$info$elapsed_s)))
metaf <- list.files(GEOM, pattern = "_meta\\.json$", full.names = TRUE)
M <- rbind_fill(lapply(metaf, function(f) { m <- fromJSON(f); data.frame(key = m$key, n = m$n, max_abs_over_g = m$max_abs_over_g, has_pts = m$has_pts,
  n_invalid_input = if (is.null(m$n_invalid_input)) 0L else m$n_invalid_input, changed = if (is.null(m$changed)) NA else m$changed,
  archetype = if (is.null(m$archetype)) NA else m$archetype, c0_c1 = if (is.null(m$c0)) NA else m$c0$c1, c0_ov = if (is.null(m$c0)) NA else m$c0$c2_overlap,
  c0_sd = if (is.null(m$c0)) NA else m$c0$c2_symdiff, roi = if (is.null(m$roi)) NA else m$roi, model = if (is.null(m$model)) NA else m$model, rep = if (is.null(m$rep)) NA else m$rep) }))
CV <- rbind_fill(lapply(list.files(COV, pattern = "\\.json$", full.names = TRUE), function(f) { j <- fromJSON(f)
  data.frame(key = j$key, v0_valid = if (is.null(j$V0)) NA else j$V0$valid, v0_ncells = if (is.null(j$V0)) NA else j$V0$n_cells_invalid_edges,
             f4_valid = if (is.null(j$F4)) NA else j$F4$valid, f4_ncells = if (is.null(j$F4)) NA else j$F4$n_cells_invalid_edges) }))
# S (perf): confronti salvati da tools/S1.3d_perf.R (rep 1)
SP <- list.files(RES, pattern = "^S_cperf_.*_(O1|O2s2|O2s8)\\.rds$", full.names = TRUE)
for (f in SP) { o <- readRDS(f)
  if (o$op == "O1") { d <- o$cmp$summary; d$key <- "S_cperf"; d$applicable <- TRUE; O1 <- rbind_fill(list(O1, d)) }
  else { d <- o$cmp$summary; d$key <- "S_cperf"; d$op <- "O2"; d$t_o2 <- o$t_op; O2 <- rbind_fill(list(O2, d))
         if (!is.null(o$cmp$o3)) { d3 <- o$cmp$o3$summary; d3$key <- "S_cperf"; d3$op <- "O3"; O2 <- rbind_fill(list(O2, d3)) } } }
for (f in SP) { o <- readRDS(f); if (o$cand == "F4" && !is.null(o$f4info$c5_n_over))
  C5 <- rbind_fill(list(C5, data.frame(key = "S_cperf", n_over = o$f4info$c5_n_over, max_ratio = o$f4info$c5_max_ratio, same_lost = o$f4info$same_lost, same_gained = o$f4info$same_gained,
                                       n_tie = o$f4info$n_tie, min_margin = o$f4info$min_margin, n_ambiguous = o$f4info$n_ambiguous, n_invalid = o$f4info$n_invalid, t_build = o$f4info$elapsed_s))) }
C5 <- C5[!duplicated(C5$key), ]
if (is.null(CV)) CV <- data.frame(key = character(0), v0_valid = logical(0), v0_ncells = numeric(0), f4_valid = logical(0), f4_ncells = numeric(0))
O1$set <- set_of(O1$key); O2$set <- set_of(O2$key)
O1 <- merge(O1, M[, c("key", "n_invalid_input", "archetype")], by = "key", all.x = TRUE); O2 <- merge(O2, M[, c("key", "n_invalid_input", "archetype")], by = "key", all.x = TRUE)
O1 <- merge(O1, CV[, c("key", "f4_valid", "v0_valid")], by = "key", all.x = TRUE)
# F4 in O1 = validatore + CoverageUnion: copertura F4 invalida (o non calcolata) = FAIL
f4 <- O1$cand == "F4" & O1$applicable %in% TRUE
O1$cov_fail <- f4 & !(O1$f4_valid %in% TRUE)
O1$fail_eff <- (O1$fail %in% TRUE) | O1$cov_fail
O1$valid_input <- O1$n_invalid_input %in% 0 | is.na(O1$n_invalid_input); O2$valid_input <- O2$n_invalid_input %in% 0 | is.na(O2$n_invalid_input)
write.csv(O1, file.path(OUTR, "S1.3d_o1_long.csv"), row.names = FALSE); write.csv(O2, file.path(OUTR, "S1.3d_o2_long.csv"), row.names = FALSE)

# ---- prestazioni (S) -------------------------------------------------------------------------------------------------
pl <- list.files(file.path(OUTD, "logs"), pattern = "^perf_.*\\.log$", full.names = TRUE)
PF <- rbind_fill(lapply(pl, function(f) { x <- readLines(f, warn = FALSE); p <- grep("^PERF,", x, value = TRUE); if (!length(p)) return(NULL)
  v <- strsplit(p[1], ",")[[1]]; rss <- as.numeric(sub(".*: ", "", grep("Maximum resident set size", x, value = TRUE)))
  data.frame(cand = v[2], op = v[3], rep = as.integer(v[4]), t_prepare = as.numeric(v[5]), t_op = as.numeric(v[6]), rss_gb = rss / 1024^2) }))
if (!is.null(PF)) { PF$t_total <- PF$t_prepare + PF$t_op
  PS <- aggregate(cbind(t_total, t_prepare, t_op, rss_gb) ~ cand + op, PF, median); PS$n_rep <- aggregate(rep ~ cand + op, PF, length)$rep
  write.csv(PS, file.path(OUTR, "S1.3d_perf.csv"), row.names = FALSE) } else PS <- NULL
tS <- function(cd, op) { if (is.null(PS)) return(NA_real_); oo <- if (op == "O1") "O1" else "O2s2"; v <- PS$t_total[PS$cand == cd & PS$op == oo]; if (length(v)) v else NA_real_ }
mS <- function(cd, op) { if (is.null(PS)) return(NA_real_); oo <- if (op == "O1") "O1" else "O2s2"; v <- PS$rss_gb[PS$cand == cd & PS$op == oo]; if (length(v)) v else NA_real_ }

# ---- controlli globali -------------------------------------------------------------------------------------------------
C6_ok <- all(M$max_abs_over_g < 1e12)
mut <- if (file.exists(file.path(RES, "mutants.rds"))) readRDS(file.path(RES, "mutants.rds")) else NULL
chf <- file.path(OUTD, "checks", c("arbiter_check_N.json", "arbiter_check_U.json")); chf <- chf[file.exists(chf)]
ach <- if (length(chf)) { tt <- lapply(chf, function(f) fromJSON(f)$totals); Reduce(function(a, b) mapply(`+`, a, b, SIMPLIFY = FALSE), tt) } else NULL
mut_ok <- !is.null(mut) && all(mut$detected) && !is.null(ach) && ach$cp4_mutant[1] == ach$cp4_mutant[2]

# ---- uscite geometriche invalide (clausola «0 invalide non riparate» della regola 1; correzione D-1, vedi dev_notes) -------
VF <- if (file.exists(file.path(OUTD, "checks", "validity_F1F3.csv"))) read.csv(file.path(OUTD, "checks", "validity_F1F3.csv")) else NULL
O1$n_invalid_out <- NA_real_
if (!is.null(VF)) { k <- O1$cand %in% c("F1", "F3"); O1$n_invalid_out[k] <- VF$n_invalid_F1F3[match(O1$key[k], VF$key)] }
k <- O1$cand == "F2"; O1$n_invalid_out[k] <- O1$still_invalid[k]
f2s <- file.path(RES, "S_cperf_F2_validity.rds"); if (file.exists(f2s)) { v <- readRDS(f2s); k <- O1$cand == "F2" & O1$key == "S_cperf"; O1$n_invalid_out[k] <- v$still_invalid; O1$n_invalid_rounded[k] <- v$n_invalid_rounded }
k <- O1$cand == "F4" & O1$applicable %in% TRUE; O1$n_invalid_out[k] <- C5$n_invalid[match(O1$key[k], C5$key)]
O1$n_invalid_out[O1$cand %in% c("V0", "F5")] <- 0          # V0: uscita = ingresso (ingressi invalidi esclusi); F5: nessuna geometria sf
inv_by <- function(cd, sets) { d <- O1[O1$cand == cd & O1$set %in% sets & O1$valid_input & O1$applicable %in% TRUE, ]; c(n = sum(d$n_invalid_out, na.rm = TRUE), na = sum(is.na(d$n_invalid_out))) }

# ---- eleggibilità ----------------------------------------------------------------------------------------------------
ELIG_SETS <- c("N", "U1", "U2", "U3", "U4", "U5", "S")       # K è dentro N; H solo dopo la decisione
elig_row <- function(op, cd) {
  if (op == "O1") { d <- O1[O1$cand == cd & O1$applicable %in% TRUE, ]
    main <- d[d$set %in% ELIG_SETS & d$valid_input, ]; cp <- d[d$set %in% c("CP1T", "CP1P") & d$valid_input, ]
    n_silent <- sum(main$silent %in% TRUE); n_exc <- sum(main$n_exceptions > 0, na.rm = TRUE) + sum(is.na(main$n_exceptions) & main$fail %in% TRUE)
    n_cov <- sum(main$cov_fail); iv <- inv_by(cd, c(ELIG_SETS, "CP1T", "CP1P")); n_inv <- iv[["n"]] + iv[["na"]]   # celle invalide (NA = non calcolato: conta)
    n_cp <- sum(cp$fail_eff); n_main <- nrow(main); hasS <- any(main$set == "S")
    ok_main <- sum(main$fail_eff) == 0
  } else { d <- O2[O2$cand == cd & O2$op == op, ]; iv <- if (cd %in% c("F1", "F2", "F3", "F4")) inv_by(cd, c(ELIG_SETS, "CP1T", "CP1P")) else c(n = 0, na = 0)
    main <- d[d$set %in% ELIG_SETS & d$valid_input, ]; cp <- d[d$set %in% c("CP1T", "CP1P") & d$valid_input, ]
    n_silent <- sum(main$silent %in% TRUE); n_exc <- sum(main$exception %in% TRUE); n_cov <- 0; n_inv <- iv[["n"]] + iv[["na"]]; n_cp <- sum(cp$fail %in% TRUE)
    n_main <- nrow(main); hasS <- any(main$set == "S"); ok_main <- sum(main$fail %in% TRUE) == 0 }
  data.frame(op = op, cand = cd, n_cases = n_main, has_S = hasS, n_silent = n_silent, n_exceptions = n_exc, n_cov_fail = n_cov, n_still_invalid = n_inv,
             n_cp1 = nrow(cp), n_cp1_fail = n_cp, mutants_ok = mut_ok,
             eligible = ok_main && n_silent == 0 && n_exc == 0 && n_cov == 0 && n_inv == 0 && n_cp == 0 && mut_ok && n_main > 0 && hasS,
             guaranteed = cd %in% GUAR[[op]] && C6_ok, t_S = tS(cd, op), rss_S = mS(cd, op))
}
EL <- do.call(rbind, lapply(names(CANDS), function(op) do.call(rbind, lapply(CANDS[[op]], function(cd) elig_row(op, cd)))))
write.csv(EL, file.path(OUTR, "S1.3d_eligibility.csv"), row.names = FALSE)

# ---- regola di decisione ---------------------------------------------------------------------------------------------
top_class <- function(op) { e <- EL[EL$op == op & EL$eligible, ]; g <- e[e$guaranteed, ]; if (nrow(g)) list(cl = "garantiti", d = g) else list(cl = "non garantiti", d = e) }
choose_op <- function(op) { tc <- top_class(op); d <- tc$d
  if (!nrow(d)) return(data.frame(op = op, class = "nessuno", chosen = NA, reason = if (op == "O1") "P0 (predicati + copertura per punti, BL-070)" else "al design dello Step 2"))
  if ("F1" %in% d$cand) { f1 <- d[d$cand == "F1", ]; ch <- d[d$cand != "F1" & d$t_S <= f1$t_S / 2 & d$rss_S <= 1.5 * f1$rss_S, ]
    if (nrow(ch)) { w <- ch[which.min(ch$t_S), ]; return(data.frame(op = op, class = tc$cl, chosen = w$cand, reason = sprintf("%.2fx più veloce di F1 su S, memoria %.2fx", f1$t_S / w$t_S, w$rss_S / f1$rss_S))) }
    return(data.frame(op = op, class = tc$cl, chosen = "F1", reason = "F1 resta (nessun candidato della classe >= 2x più veloce con memoria <= 1.5x)")) }
  w <- d[which.min(d$t_S), ]; data.frame(op = op, class = tc$cl, chosen = w$cand, reason = "F1 fuori dalla classe più alta: il più veloce") }
DEC <- do.call(rbind, lapply(c("O1", "O2", "O3"), choose_op))
# un solo helper (BL-081)
single <- NA; Tcomb <- NA
if (all(!is.na(DEC$chosen)) && length(unique(DEC$chosen)) > 1) {
  adm <- Reduce(intersect, lapply(c("O1", "O2", "O3"), function(op) top_class(op)$d$cand))
  Tcomb <- tS(DEC$chosen[DEC$op == "O1"], "O1") + tS(DEC$chosen[DEC$op == "O2"], "O2")
  if (length(adm)) { Tc <- vapply(adm, function(cd) tS(cd, "O1") + tS(cd, "O2"), 0); best <- adm[which.min(Tc)]
    if (!(Tcomb <= min(Tc) / 2)) { single <- best; DEC$reason <- paste0(DEC$reason, sprintf(" | un solo helper: %s (T %.1f s vs combinazione %.1f s)", best, min(Tc), Tcomb)); DEC$chosen_final <- best } } }
if (is.null(DEC$chosen_final)) DEC$chosen_final <- DEC$chosen
write.csv(DEC, file.path(OUTR, "S1.3d_decision.csv"), row.names = FALSE)
print(EL, row.names = FALSE); print(DEC, row.names = FALSE)

# ---- verdetti (C, previsioni, CP) ------------------------------------------------------------------------------------
V <- list(); vr <- function(section, id, cand, metric, value, threshold, outcome, prediction = NA, confirmed = NA)
  V[[length(V) + 1]] <<- data.frame(section = section, id = id, cand = cand, metric = metric, value = as.character(value), threshold = as.character(threshold),
                                    outcome = outcome, prediction = as.character(prediction), confirmed = confirmed)
pf <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
# C-0
ns13 <- read.csv("results/S1.3/S1.3_null_summary.csv"); ns13$key <- sprintf("%s_%s_%s_%02d", ns13$archetype, ns13$roi_id, ns13$model, ns13$rep)
MN <- merge(M[set_of(M$key) == "N", ], ns13[, c("key", "n", "c1", "c2_overlap", "c2_symdiff")], by = "key", suffixes = c("", ".s13"))
rel <- function(x, y) abs(x - y) / pmax(abs(y), 1e-300)
c0ok <- MN$n == MN$n.s13 & (rel(MN$c0_c1, MN$c1) <= 1e-12 | MN$c0_c1 == MN$c1) & (rel(MN$c0_ov, MN$c2_overlap) <= 1e-12 | MN$c0_ov == MN$c2_overlap) & (rel(MN$c0_sd, MN$c2_symdiff) <= 1e-12 | MN$c0_sd == MN$c2_symdiff)
vr("C", "C-0", "-", "nulli riprodotti (n, c1, c2 rel <= 1e-12)", sprintf("%d/%d", sum(c0ok), nrow(MN)), "1800/1800", pf(sum(c0ok) == 1800 && nrow(MN) == 1800), "1800/1800", sum(c0ok) == 1800)
tst <- read.csv("results/S1.3d/S1.3d_test.csv"); vr("C", "C-1", "-", "test unitari", sprintf("%d/%d", sum(tst$ok), nrow(tst)), "100 %", pf(all(tst$ok)), "100 %", all(tst$ok))
if (!is.null(ach)) for (t in c("c2a", "c2b", "c2c")) vr("C", paste0("C-2", substr(t, 3, 3)), "arbitro", c(c2a = "cella ∩ regione: indice vs tutti i lati", c2b = "coppie: Green vs S-H", c2c = "pezzi O2: S-H vs Green")[[t]],
  sprintf("%d/%d", ach[[t]][1], ach[[t]][2]), "100 %", pf(ach[[t]][1] == ach[[t]][2] && ach[[t]][2] > 0), "100 %", ach[[t]][1] == ach[[t]][2])
if (!is.null(mut)) for (q in seq_len(nrow(mut))) vr("C", "C-3", mut$mutant[q], mut$detail[q], mut$detected[q], "rilevato", pf(mut$detected[q]), "rilevato", mut$detected[q])
if (!is.null(ach)) vr("C", "C-3", "M-arb", "vertice tolto -> area esatta cambia", sprintf("%d/%d", ach$cp4_mutant[1], ach$cp4_mutant[2]), "100 %", pf(ach$cp4_mutant[1] == ach$cp4_mutant[2]), "100 %", ach$cp4_mutant[1] == ach$cp4_mutant[2])
ch <- M$changed[set_of(M$key) %in% c("U1", "U2", "U4", "CP1T", "CP1P") | grepl("^U3_ulp", M$key)]
vr("C", "C-4", "-", "perturbazioni che cambiano l'input", sprintf("%d/%d", sum(ch %in% TRUE), length(ch)), "100 %", pf(all(ch %in% TRUE)), "100 %", all(ch %in% TRUE))
c5 <- C5[!is.na(C5$n_over), ]; c5ok <- c5$n_over == 0 & c5$same_lost & c5$same_gained
vr("C", "C-5", "F4", "F4 vs 04b3: |Δa| <= tau_i e stessi frammenti (casi; parità)", sprintf("%d/%d; parità %d; errori di costruzione %d", sum(c5ok), nrow(c5), sum(c5$n_tie), sum(!is.na(C5$error))),
   "100 % fuori dalle parità", pf(all(c5ok | c5$n_tie > 0) && all(is.na(C5$error))), "100 %", all(c5ok))
vr("C", "C-6", "-", "max |coord| / g", sprintf("%.3e", max(M$max_abs_over_g)), "< 1e12", pf(C6_ok), "< 1e12", C6_ok)
# previsioni
o1N <- O1[O1$set == "N", ]; Kk <- NULL
pred <- function(id, cand, metric, value, ok, pr) vr("P", id, cand, metric, value, pr, if (isTRUE(ok)) "confermata" else "smentita", pr, isTRUE(ok))
v0 <- o1N[o1N$cand == "V0", ]
vr("P", "P-V0a", "V0", "nulli N con errore silenzioso O1 (unione / coppie)", sprintf("%d (unione %d, coppie %d)", sum(v0$silent %in% TRUE), sum(v0$union_viol %in% TRUE), sum(v0$pair_viol_any > 0, na.rm = TRUE)),
   ">= 15 unioni e >= 11 coppie", if (sum(v0$union_viol %in% TRUE) >= 15 && sum(v0$pair_viol_any, na.rm = TRUE) >= 11) "confermata" else "smentita", "15 unioni, >= 11 coppie",
   sum(v0$union_viol %in% TRUE) >= 15 && sum(v0$pair_viol_any, na.rm = TRUE) >= 11)
v02 <- O2[O2$set == "N" & O2$cand == "V0", ]; vr("P", "P-V0b", "V0", "violazioni O2/O3 su N", sum(v02$n_viol, na.rm = TRUE), "0", if (sum(v02$n_viol, na.rm = TRUE) == 0) "confermata" else "smentita", "0", sum(v02$n_viol, na.rm = TRUE) == 0)
for (cd in c("F1", "F2", "F3", "F4", "F5")) { e <- EL[EL$cand == cd & EL$op == "O1", ]
  vr("P", paste0("P-", cd), cd, "eleggibile O1 (errori silenziosi, eccezioni, CP-1)", sprintf("%s (silenziosi %d, eccezioni %d, CP-1 %d)", e$eligible, e$n_silent, e$n_exceptions, e$n_cp1_fail), "eleggibile",
     if (e$eligible) "confermata" else "smentita", "eleggibile", e$eligible) }
for (cd in c("F1", "F6")) for (op in c("O2", "O3")) { e <- EL[EL$cand == cd & EL$op == op, ]
  vr("P", paste0("P-", cd, "-", op), cd, paste("eleggibile", op), sprintf("%s (silenziosi %d)", e$eligible, e$n_silent), "eleggibile", if (e$eligible) "confermata" else "smentita", "eleggibile", e$eligible) }
f2N <- O1[O1$cand == "F2" & O1$set == "N", ]; vr("P", "P-F2inv", "F2", "celle invalide dopo l'arrotondamento in N (nulli; ancora invalide)", sprintf("%d; %d", sum(f2N$n_invalid_rounded, na.rm = TRUE), sum(f2N$still_invalid, na.rm = TRUE)),
   "2 rilevate e riparate", if (sum(f2N$n_invalid_rounded, na.rm = TRUE) == 2 && sum(f2N$still_invalid, na.rm = TRUE) == 0) "confermata" else "smentita", "2; 0",
   sum(f2N$n_invalid_rounded, na.rm = TRUE) == 2 && sum(f2N$still_invalid, na.rm = TRUE) == 0)
cvN <- CV[set_of(CV$key) == "N", ]
vr("P", "P-F4cov", "F4", "copertura F4 valida (N)", sprintf("%d/%d", sum(cvN$f4_valid %in% TRUE), nrow(cvN)), "100 %", if (all(cvN$f4_valid %in% TRUE)) "confermata" else "smentita", "100 %", all(cvN$f4_valid %in% TRUE))
vr("P", "P-V0cov", "V0", "copertura V0 invalida (N)", sprintf("%d/%d (%.1f %%)", sum(cvN$v0_valid %in% FALSE), nrow(cvN), 100 * mean(cvN$v0_valid %in% FALSE)), ">= 80 %",
   if (mean(cvN$v0_valid %in% FALSE) >= 0.8) "confermata" else "smentita", ">= 80 %", mean(cvN$v0_valid %in% FALSE) >= 0.8)
# CP-1, CP-2
for (cd in CANDS$O1) { d <- O1[O1$cand == cd & O1$set %in% c("CP1T", "CP1P") & O1$applicable %in% TRUE, ]
  vr("CP", "CP-1", cd, "FAIL su input validi (O1); input invalidi esclusi", sprintf("%d/%d; esclusi %d", sum(d$fail_eff[d$valid_input]), sum(d$valid_input), sum(!d$valid_input)), "0",
     if (sum(d$fail_eff[d$valid_input]) == 0) "PASS" else "FAIL") }
u4 <- O1[O1$set == "U4" & O1$applicable %in% TRUE, ]
for (cd in unique(u4$cand)) { d <- u4[u4$cand == cd, ]; vr("CP", "CP-2", cd, "U4 a cavallo della griglia: casi oltre tau (max errore/tau coppie)", sprintf("%d/%d (%.3g)", sum(d$fail_eff), nrow(d), max(d$pair_max_ratio, na.rm = TRUE)), "0",
     if (sum(d$fail_eff) == 0) "PASS" else "FAIL") }
if (!is.null(ach)) vr("CP", "CP-4", "arbitro", "placebo (ordine permutato: input cambiato, verità identica)", sprintf("%d/%d", ach$cp4_placebo[1], ach$cp4_placebo[2]), "100 %", pf(ach$cp4_placebo[1] == ach$cp4_placebo[2]))
# CP-3 (holdout, solo dopo la decisione)
if (HOLD) for (op in c("O1", "O2", "O3")) { cd <- DEC$chosen_final[DEC$op == op]; if (is.na(cd)) next
  if (op == "O1") { d <- O1[O1$set == "H" & O1$cand == cd, ]; nf <- sum(d$fail_eff) } else { d <- O2[O2$set == "H" & O2$cand == cd & O2$op == op, ]; nf <- sum(d$fail %in% TRUE) }
  vr("CP", "CP-3", cd, paste("H (300 nulli nuovi):", op, "casi con errore"), sprintf("%d/%d", nf, nrow(d)), "0", if (nf == 0 && nrow(d) > 0) "PASS" else "FAIL") }
V <- do.call(rbind, V); write.csv(V, file.path(OUTR, "S1.3d_verdicts.csv"), row.names = FALSE); print(V[, c("section", "id", "cand", "value", "outcome")], row.names = FALSE)
