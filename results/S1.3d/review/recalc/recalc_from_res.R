# results/S1.3d/review/recalc/recalc_from_res.R — revisione avversariale «ricalcolo» di S1.3d.
# Ricalcolo INDIPENDENTE (non usa tools/S1.3d_analysis.R né i CSV aggregati come input) dai grezzi:
#   res/*_o1.rds, *_o12.rds, S_cperf_*.rds, coverage/*.json, geom/*_meta.json, checks/*, logs/perf_*.log, b1/verdicts_*.csv
# I CSV ufficiali (eligibility, decision, verdicts, b1) sono letti SOLO alla fine, per il confronto.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.3d/review/recalc/recalc_from_res.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(jsonlite); library(parallel) })
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; HEAVY <- file.path(D, "review_recalc"); REV <- "results/S1.3d/review/recalc"
g <- 1e-6; NCR <- as.integer(Sys.getenv("RC_CORES", "1"))
rk <- function(f) readLines(file.path(D, f))
KN <- rk("keys_N.txt"); KH <- rk("keys_H.txt"); KU <- rk("keys_U.txt"); KCP <- rk("keys_CP1.txt"); KK <- rk("keys_K.txt"); KO2 <- rk("keys_O2.txt"); KUo2 <- rk("keys_U_o2.txt")
cat("== insiemi: N", length(KN), "uniq", length(unique(KN)), "| H", length(KH), "| U", length(KU), "| CP1", length(KCP), "| K", length(KK), "| O2", length(KO2), "| U_o2", length(KUo2), "\n")
stopifnot(!anyDuplicated(c(KN, KH, KU, KCP)))
cat("K in N:", all(KK %in% KN), " O2 in N:", all(KO2 %in% KN), " K in O2:", all(KK %in% KO2), "\n")
repN <- as.integer(sub(".*_", "", KN)); repH <- as.integer(sub(".*_", "", KH))
cat("rep N:", paste(range(repN), collapse = ".."), " rep H:", paste(sort(unique(repH)), collapse = ","), " H per modello:", paste(names(table(sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", KH))), table(sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", KH)), collapse = " "), "\n")
subset_of <- function(k) { s <- rep(NA_character_, length(k))
  s[k %in% KN] <- "N"; s[k %in% KH] <- "H"; s[k == "S_cperf"] <- "S"
  for (u in paste0("U", 1:5)) s[k %in% KU & startsWith(k, paste0(u, "_"))] <- u
  s[k %in% KCP & startsWith(k, "CP1T_")] <- "CP1T"; s[k %in% KCP & startsWith(k, "CP1P_")] <- "CP1P"; s }
ALLK <- c(KN, KH, KU, KCP)
cat("chiavi senza sottoinsieme:", sum(is.na(subset_of(ALLK))), "\n")

# ---- meta ------------------------------------------------------------------------------------------------------------
meta <- mclapply(c(ALLK, "S_cperf"), function(k) fromJSON(file.path(D, "geom", paste0(k, "_meta.json"))), mc.cores = NCR)
names(meta) <- c(ALLK, "S_cperf")
M <- data.frame(key = names(meta), n = vapply(meta, function(m) as.numeric(m$n), 0), Preg = vapply(meta, function(m) m$Preg, 0), sumP = vapply(meta, function(m) sum(m$P), 0),
                lenP = vapply(meta, function(m) length(m$P), 0), maxg = vapply(meta, function(m) m$max_abs_over_g, 0), has_pts = vapply(meta, function(m) isTRUE(m$has_pts), TRUE),
                n_inv_in = vapply(meta, function(m) if (is.null(m$n_invalid_input)) NA_real_ else as.numeric(m$n_invalid_input), 0),
                changed = vapply(meta, function(m) if (is.null(m$changed)) NA else as.logical(m$changed), NA), stringsAsFactors = FALSE)
M$set <- subset_of(M$key)
cat("meta: lenP == n:", all(M$lenP == M$n), "; max |coord|/g =", format(max(M$maxg), digits = 4), "(C-6 < 1e12:", max(M$maxg) < 1e12, ")\n")
cat("n_invalid_input disponibile per insieme:\n"); print(table(M$set, is.na(M$n_inv_in)))
# input valido: per gli insiemi costruiti n_invalid_input == 0; per N/H/S non registrato -> validità ricalcolata in recalc_independent.R (input_invalid)
IV <- if (file.exists(file.path(HEAVY, "input_validity.csv"))) read.csv(file.path(HEAVY, "input_validity.csv")) else NULL
M$n_inv_in2 <- M$n_inv_in
if (!is.null(IV)) { k <- is.na(M$n_inv_in2); M$n_inv_in2[k] <- IV$n_invalid_input[match(M$key[k], IV$key)] }
M$valid_input <- M$n_inv_in2 %in% 0
cat("input validi per insieme (NA = non verificato):\n"); print(table(M$set, ifelse(is.na(M$n_inv_in2), "NA", ifelse(M$n_inv_in2 == 0, "valido", "invalido"))))

# ---- res: O1/O2 per caso ---------------------------------------------------------------------------------------------
rf <- function(k) { a <- file.path(D, "res", paste0(k, "_o12.rds")); b <- file.path(D, "res", paste0(k, "_o1.rds")); c(o12 = file.exists(a), o1 = file.exists(b)) }
FX <- t(vapply(ALLK, rf, c(o12 = TRUE, o1 = TRUE)))
cat("res: chiavi con entrambi i file", sum(FX[, 1] & FX[, 2]), "; con nessuno", sum(!FX[, 1] & !FX[, 2]), "\n")
expect_o2 <- ALLK %in% c(KO2, KUo2, KCP, KH)
cat("O2 atteso ma _o12 assente:", sum(expect_o2 & !FX[, 1]), "; _o12 presente ma non atteso:", sum(!expect_o2 & FX[, 1]), "\n")
num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)
lg <- function(x) if (is.null(x) || !length(x)) NA else as.logical(x)
one <- function(k) {
  f <- file.path(D, "res", paste0(k, if (FX[k, "o12"]) "_o12" else "_o1", ".rds")); o <- readRDS(f); m <- meta[[k]]
  stopifnot(identical(o$key, k))
  d <- o$o1; rows1 <- lapply(seq_len(nrow(d)), function(q) { r <- d[q, , drop = FALSE]; gv <- function(nm) if (nm %in% names(r)) r[[nm]] else NA
    data.frame(key = k, cand = r$cand, applicable = lg(gv("applicable")), union_err = num(gv("union_err")), symd_err = num(gv("symd_err")), tauU_st = num(gv("tauU")),
               union_viol_st = lg(gv("union_viol")), symd_viol_st = lg(gv("symd_viol")), pv_ij = num(gv("pair_viol_ij")), pv_ji = num(gv("pair_viol_ji")), pv_any = num(gv("pair_viol_any")),
               cell_viol = num(gv("cell_viol")), n_exc = num(gv("n_exceptions")), silent_st = lg(gv("silent")), fail_st = lg(gv("fail")), err_msg = as.character(gv("error")),
               still_invalid = num(gv("still_invalid")), n_inv_rounded = num(gv("n_invalid_rounded")), pmax_ratio = num(gv("pair_max_ratio")), cmax_ratio = num(gv("cell_max_ratio")), stringsAsFactors = FALSE) })
  O1 <- do.call(rbind, rows1)
  # violazioni registrate: conteggio per tipo e tolleranza ricalcolata da meta$P
  vchk <- data.frame(key = k, cand = O1$cand, n_vij = 0, n_vji = 0, n_vcell = 0, tau_bad = 0)
  for (cd in O1$cand) { v <- o$viol[[paste0("o1_", cd)]]; if (is.null(v)) next; q <- which(vchk$cand == cd)
    vchk$n_vij[q] <- sum(v$what == "pair_ij"); vchk$n_vji[q] <- sum(v$what == "pair_ji"); vchk$n_vcell[q] <- sum(v$what == "cell")
    pp <- v$what %in% c("pair_ij", "pair_ji"); cc <- v$what == "cell"
    tre <- rep(NA_real_, nrow(v)); tre[pp] <- g * (m$P[v$i[pp]] + m$P[v$j[pp]]); tre[cc] <- g * m$P[v$i[cc]]
    vchk$tau_bad[q] <- sum(abs(tre[pp | cc] - v$tau[pp | cc]) > 1e-15 * pmax(1, abs(v$tau[pp | cc]))) }
  O2 <- NULL
  if (length(o$o2) && is.data.frame(o$o2)) { e <- o$o2
    O2 <- data.frame(key = k, cand = e$cand, op = ifelse(e$op == "O3", "O3", paste0("O2s", e$s)), n_truth = num(e$n_truth), n_cand = num(e$n_cand), n_viol = num(e$n_viol), n_missing = num(e$n_missing),
                     n_spurious = num(e$n_spurious), max_ratio = num(e$max_ratio), exception = as.logical(e$exception), fail_st = as.logical(e$fail), silent_st = as.logical(e$silent), stringsAsFactors = FALSE) }
  C5 <- if (is.null(o$c5)) NULL else if (!is.null(o$c5$error)) data.frame(key = k, error = substr(o$c5$error, 1, 100), n_over = NA, same_lost = NA, same_gained = NA, n_tie = NA, n_invalid = NA, n_ambiguous = NA) else
    data.frame(key = k, error = NA, n_over = o$c5$n_over, same_lost = o$c5$same_lost, same_gained = o$c5$same_gained, n_tie = o$c5$info$n_tie, n_invalid = o$c5$info$n_invalid, n_ambiguous = o$c5$info$n_ambiguous)
  list(O1 = O1, O2 = O2, C5 = C5, V = vchk, mt = file.mtime(f))
}
R <- mclapply(ALLK, function(k) tryCatch(one(k), error = function(e) paste("ERR", k, conditionMessage(e))), mc.cores = NCR, mc.preschedule = FALSE)
bad <- vapply(R, is.character, TRUE); if (any(bad)) { print(unlist(R[bad])); stop("errori di lettura") }
O1 <- do.call(rbind, lapply(R, `[[`, "O1")); O2 <- do.call(rbind, lapply(R, `[[`, "O2")); C5 <- do.call(rbind, lapply(R, `[[`, "C5")); VC <- do.call(rbind, lapply(R, `[[`, "V"))
MT <- data.frame(key = ALLK, mtime = as.POSIXct(vapply(R, function(r) as.numeric(r$mt), 0), origin = "1970-01-01"))

# ---- S (S_cperf_*) ----------------------------------------------------------------------------------------------------
sf_ <- list.files(file.path(D, "res"), pattern = "^S_cperf_.*\\.rds$"); cat("file S:", paste(sf_, collapse = " "), "\n")
mS <- meta[["S_cperf"]]
for (f in sf_[grepl("_(O1|O2s2|O2s8)\\.rds$", sf_)]) { o <- readRDS(file.path(D, "res", f))
  if (o$op == "O1") { r <- o$cmp$summary
    O1 <- rbind(O1, data.frame(key = "S_cperf", cand = o$cand, applicable = TRUE, union_err = r$union_err, symd_err = r$symd_err, tauU_st = r$tauU, union_viol_st = r$union_viol, symd_viol_st = r$symd_viol,
      pv_ij = r$pair_viol_ij, pv_ji = r$pair_viol_ji, pv_any = r$pair_viol_any, cell_viol = r$cell_viol, n_exc = r$n_exceptions, silent_st = r$silent, fail_st = r$fail, err_msg = NA,
      still_invalid = NA, n_inv_rounded = NA, pmax_ratio = r$pair_max_ratio, cmax_ratio = r$cell_max_ratio, stringsAsFactors = FALSE))
    if (o$cand == "F4") C5 <- rbind(C5, data.frame(key = "S_cperf", error = NA, n_over = o$f4info$c5_n_over, same_lost = o$f4info$same_lost, same_gained = o$f4info$same_gained,
                                                   n_tie = o$f4info$n_tie, n_invalid = o$f4info$n_invalid, n_ambiguous = o$f4info$n_ambiguous))
  } else { add <- function(r, op) data.frame(key = "S_cperf", cand = o$cand, op = op, n_truth = r$n_truth, n_cand = r$n_cand, n_viol = r$n_viol, n_missing = r$n_missing, n_spurious = r$n_spurious,
                                          max_ratio = r$max_ratio, exception = r$exception, fail_st = r$fail, silent_st = r$silent, stringsAsFactors = FALSE)
    O2 <- rbind(O2, add(o$cmp$summary, o$op)); if (!is.null(o$cmp$o3)) O2 <- rbind(O2, add(o$cmp$o3$summary, "O3")) }
}
f2v <- readRDS(file.path(D, "res", "S_cperf_F2_validity.rds")); O1$still_invalid[O1$key == "S_cperf" & O1$cand == "F2"] <- f2v$still_invalid
O1$n_inv_rounded[O1$key == "S_cperf" & O1$cand == "F2"] <- f2v$n_invalid_rounded
cat("S: righe O1", sum(O1$key == "S_cperf"), "(", paste(O1$cand[O1$key == "S_cperf"], collapse = ","), ") righe O2/O3", sum(O2$key == "S_cperf"), "\n")
O1$set <- subset_of(O1$key); O2$set <- subset_of(O2$key)
C5$set <- subset_of(C5$key)

# ---- ri-derivazione dei flag O1 dalle quantità grezze (tau da meta) -------------------------------------------------
O1$tauU <- g * (M$Preg[match(O1$key, M$key)] + M$sumP[match(O1$key, M$key)])
O1$tau_rel <- abs(O1$tauU - O1$tauU_st) / O1$tauU
O1$f4_build_fail <- O1$cand == "F4" & O1$applicable %in% TRUE & is.na(O1$union_err) & is.na(O1$n_exc) & O1$fail_st %in% TRUE
O1$uv <- is.na(O1$union_err) | O1$union_err > O1$tauU; O1$sv <- is.na(O1$symd_err) | O1$symd_err > O1$tauU
O1$exc <- ifelse(is.na(O1$n_exc), as.numeric(O1$fail_st %in% TRUE & O1$applicable %in% TRUE), O1$n_exc)
O1$viol_any <- O1$uv | O1$sv | (O1$pv_any %in% NA) | O1$pv_any > 0 | (O1$cell_viol %in% NA) | O1$cell_viol > 0
O1$viol_any[is.na(O1$viol_any)] <- TRUE
O1$silent <- O1$viol_any & O1$exc == 0 & O1$applicable %in% TRUE
O1$fail <- (O1$viol_any | O1$exc > 0) & O1$applicable %in% TRUE
ap <- O1$applicable %in% TRUE & !O1$f4_build_fail
cat("\n== coerenza flag O1 (righe applicabili", sum(ap), ")\n")
cat("tauU ricalcolato vs registrato: max rel", format(max(O1$tau_rel[ap], na.rm = TRUE), digits = 3), "\n")
cat("union_viol diversi:", sum(O1$uv[ap] != O1$union_viol_st[ap], na.rm = TRUE), " symd_viol diversi:", sum(O1$sv[ap] != O1$symd_viol_st[ap], na.rm = TRUE),
    " silent diversi:", sum(O1$silent[ap] != O1$silent_st[ap], na.rm = TRUE), " fail diversi:", sum(O1$fail[ap] != O1$fail_st[ap], na.rm = TRUE), "\n")
cat("pv_any vs |ij ∪ ji|: pv_any < max(ij,ji):", sum(O1$pv_any[ap] < pmax(O1$pv_ij[ap], O1$pv_ji[ap]), na.rm = TRUE), " pv_any > ij+ji:", sum(O1$pv_any[ap] > O1$pv_ij[ap] + O1$pv_ji[ap], na.rm = TRUE), "\n")
cat("max_ratio coppie > 1 con 0 violazioni:", sum(O1$pmax_ratio[ap] > 1 & O1$pv_any[ap] == 0, na.rm = TRUE), "; max_ratio celle > 1 con 0 violazioni:", sum(O1$cmax_ratio[ap] > 1 & O1$cell_viol[ap] == 0, na.rm = TRUE), "\n")
VC2 <- merge(VC, O1[, c("key", "cand", "pv_ij", "pv_ji", "cell_viol")], by = c("key", "cand"))
cat("tabelle viol: righe pair_ij != pv_ij:", sum(VC2$n_vij != VC2$pv_ij, na.rm = TRUE), " pair_ji:", sum(VC2$n_vji != VC2$pv_ji, na.rm = TRUE), " cell:", sum(VC2$n_vcell != VC2$cell_viol, na.rm = TRUE),
    " tau delle violazioni != g*(P_i+P_j) / g*P_i:", sum(VC2$tau_bad), "\n")
cat("F4: costruzione fallita (righe):", sum(O1$f4_build_fail), " errori c5:", sum(!is.na(C5$error)), "\n")
cat("righe O1 per cand x applicabile:\n"); print(table(O1$cand, O1$applicable, useNA = "ifany"))

# ---- O2/O3 ------------------------------------------------------------------------------------------------------------
O2$fail <- O2$exception %in% TRUE | is.na(O2$n_viol) | O2$n_viol > 0
cat("\n== coerenza flag O2/O3: fail diversi", sum(O2$fail != O2$fail_st), "; silent != fail (senza eccezione):", sum(O2$silent_st != (O2$fail & !O2$exception %in% TRUE)),
    "; max_ratio > 1 con 0 violazioni:", sum(O2$max_ratio > 1 & O2$n_viol == 0, na.rm = TRUE), "; n_cand NA:", sum(is.na(O2$n_cand)), "\n")
exp_c <- function(k) { h <- M$has_pts[M$key == k]; c("V0", "F1", "F2", "F3", if (h) "F4", "F5", "F6") }
o2k <- unique(O2$key); miss <- 0
for (k in o2k) for (cd in exp_c(k)) for (op in c("O2s2", "O2s8", "O3")) if (!any(O2$key == k & O2$cand == cd & O2$op == op)) miss <- miss + 1
cat("righe O2/O3 attese ma assenti (F4 solo con generatori):", miss, " su", length(o2k), "casi\n")
print(table(O2$set, O2$op))

# ---- copertura -----------------------------------------------------------------------------------------------------------
CV <- do.call(rbind, mclapply(c(ALLK, "S_cperf"), function(k) { f <- file.path(D, "coverage", paste0(k, ".json")); if (!file.exists(f)) return(data.frame(key = k, v0 = NA, f4 = NA, v0n = NA, present = FALSE))
  j <- fromJSON(f); data.frame(key = k, v0 = if (is.null(j$V0)) NA else j$V0$valid, f4 = if (is.null(j$F4)) NA else j$F4$valid, v0n = if (is.null(j$V0)) NA else j$V0$n_cells_invalid_edges, present = TRUE) }, mc.cores = NCR))
CV$set <- subset_of(CV$key); cat("\ncopertura: json mancanti", sum(!CV$present), "\n")
O1$f4_cov_valid <- CV$f4[match(O1$key, CV$key)]
O1$cov_fail <- O1$cand == "F4" & O1$applicable %in% TRUE & !(O1$f4_cov_valid %in% TRUE)
O1$fail_eff <- O1$fail | O1$cov_fail

# ---- celle invalide (D-1) --------------------------------------------------------------------------------------------
VF <- read.csv(file.path(D, "checks", "validity_F1F3.csv"))
cat("validity_F1F3.csv: chiavi", nrow(VF), "; chiavi H presenti:", sum(VF$key %in% KH), "; S:", "S_cperf" %in% VF$key, "\n")
O1$n_inv_out <- NA_real_
k <- O1$cand %in% c("F1", "F3"); O1$n_inv_out[k] <- VF$n_invalid_F1F3[match(O1$key[k], VF$key)]
k <- O1$cand == "F2"; O1$n_inv_out[k] <- O1$still_invalid[k]
k <- O1$cand == "F4" & O1$applicable %in% TRUE; O1$n_inv_out[k] <- C5$n_invalid[match(O1$key[k], C5$key)]
k <- O1$cand %in% c("V0", "F5"); O1$n_inv_out[k] <- 0
O1$valid_input <- M$valid_input[match(O1$key, M$key)]; O2$valid_input <- M$valid_input[match(O2$key, M$key)]
saveRDS(list(O1 = O1, O2 = O2, C5 = C5, CV = CV, M = M, MT = MT), file.path(HEAVY, "recalc_tables.rds"))

# ---- perf ---------------------------------------------------------------------------------------------------------------
pl <- list.files(file.path(D, "logs"), pattern = "^perf_(V0|F[1-6])_(O1|O2s2|O2s8)_[0-9]+\\.log$", full.names = TRUE)
PF <- do.call(rbind, lapply(pl, function(f) { x <- readLines(f, warn = FALSE); b <- strsplit(sub("\\.log$", "", sub("^perf_", "", basename(f))), "_")[[1]]
  p <- grep("^PERF,", x, value = TRUE); rss <- grep("Maximum resident set size", x, value = TRUE); wall <- grep("Elapsed \\(wall clock\\)", x, value = TRUE)
  v <- if (length(p)) strsplit(p[1], ",")[[1]] else rep(NA, 6)
  data.frame(cand = b[1], op = b[2], rep = as.integer(b[3]), has_perf = length(p) > 0, perf_cand = v[2], perf_op = v[3], perf_rep = v[4], t_prep = as.numeric(v[5]), t_op = as.numeric(v[6]),
             rss_gb = if (length(rss)) as.numeric(sub(".*: ", "", rss[1])) / 1024^2 else NA, wall = if (length(wall)) sub(".*: ", "", wall[1]) else NA, size = file.size(f), mtime = format(file.mtime(f)), stringsAsFactors = FALSE) }))
cat("\n== perf: log", nrow(PF), "; con riga PERF", sum(PF$has_perf), "; PERF incoerente col nome file:", sum(PF$has_perf & (PF$perf_cand != PF$cand | PF$perf_op != PF$op | as.integer(PF$perf_rep) != PF$rep)), "\n")
PF$t_tot <- PF$t_prep + PF$t_op
PS <- do.call(rbind, lapply(split(PF[PF$has_perf, ], paste(PF$cand[PF$has_perf], PF$op[PF$has_perf])), function(d)
  data.frame(cand = d$cand[1], op = d$op[1], n_rep = nrow(d), t_med = median(d$t_tot), t_min = min(d$t_tot), t_max = max(d$t_tot), rss_med = median(d$rss_gb), rss_max = max(d$rss_gb))))
print(PS, row.names = FALSE); print(PF[!PF$has_perf, c("cand", "op", "rep", "size", "mtime", "wall")], row.names = FALSE)
write.csv(PF, file.path(REV, "recalc_perf_logs.csv"), row.names = FALSE); write.csv(PS, file.path(REV, "recalc_perf_median.csv"), row.names = FALSE)
tS <- function(cd, op) { oo <- if (op == "O1") "O1" else "O2s2"; v <- PS$t_med[PS$cand == cd & PS$op == oo]; if (length(v)) v else NA_real_ }
mS <- function(cd, op) { oo <- if (op == "O1") "O1" else "O2s2"; v <- PS$rss_med[PS$cand == cd & PS$op == oo]; if (length(v)) v else NA_real_ }

# ---- mutanti, C-6 ------------------------------------------------------------------------------------------------------
mut <- readRDS(file.path(D, "res", "mutants.rds")); ach <- lapply(c("N", "U"), function(s) fromJSON(file.path(D, "checks", sprintf("arbiter_check_%s.json", s)))$totals)
marb <- Reduce(`+`, lapply(ach, function(a) a$cp4_mutant)); plac <- Reduce(`+`, lapply(ach, function(a) a$cp4_placebo))
c2 <- sapply(c("c2a", "c2b", "c2c"), function(t) Reduce(`+`, lapply(ach, function(a) a[[t]])))
mut_ok <- all(mut$detected) && marb[1] == marb[2]
cat("\nmutanti C-3:", paste(mut$mutant, mut$detected, collapse = "; "), "| M-arb", marb[1], "/", marb[2], "| CP-4 placebo", plac[1], "/", plac[2], "| C-2", paste(colnames(c2), c2[1, ], "/", c2[2, ], collapse = "; "), "\n")
C6 <- max(M$maxg) < 1e12

# ---- eleggibilità (regola 1 della pre-registrazione) -----------------------------------------------------------------
MAIN <- c("N", "U1", "U2", "U3", "U4", "U5", "S"); CPS <- c("CP1T", "CP1P")
GUAR <- list(O1 = c("F1", "F5", "F4"), O2 = c("F1", "F5"), O3 = c("F1", "F5"))
CANDS <- list(O1 = c("V0", "F1", "F2", "F3", "F4", "F5"), O2 = c("V0", "F1", "F2", "F3", "F4", "F5", "F6"), O3 = c("V0", "F1", "F2", "F3", "F4", "F5", "F6"))
invc <- function(cd, sets) { d <- O1[O1$cand == cd & O1$set %in% sets & O1$valid_input %in% TRUE & O1$applicable %in% TRUE, ]; c(n = sum(d$n_inv_out, na.rm = TRUE), na = sum(is.na(d$n_inv_out)), ncase = sum(d$n_inv_out > 0, na.rm = TRUE)) }
EL <- list()
for (op in names(CANDS)) for (cd in CANDS[[op]]) {
  if (op == "O1") { d <- O1[O1$cand == cd & O1$applicable %in% TRUE, ]; mn <- d[d$set %in% MAIN & d$valid_input %in% TRUE, ]; cp <- d[d$set %in% CPS & d$valid_input %in% TRUE, ]
    n_sil <- sum(mn$silent); n_exc <- sum(mn$exc > 0); n_cov <- sum(mn$cov_fail); n_fail <- sum(mn$fail_eff); n_cpf <- sum(cp$fail_eff)
  } else { d <- O2[O2$cand == cd & O2$op == (if (op == "O2") c("O2s2", "O2s8") else "O3")[1], ]
    if (op == "O2") d <- O2[O2$cand == cd & O2$op %in% c("O2s2", "O2s8"), ]
    mn <- d[d$set %in% MAIN & d$valid_input %in% TRUE, ]; cp <- d[d$set %in% CPS & d$valid_input %in% TRUE, ]
    n_sil <- sum(mn$fail & !mn$exception %in% TRUE); n_exc <- sum(mn$exception %in% TRUE); n_cov <- 0; n_fail <- sum(mn$fail); n_cpf <- sum(cp$fail) }
  iv_main <- if (cd %in% c("F1", "F2", "F3", "F4")) invc(cd, MAIN) else c(n = 0, na = 0, ncase = 0)
  iv_cp <- if (cd %in% c("F1", "F2", "F3", "F4")) invc(cd, CPS) else c(n = 0, na = 0, ncase = 0)
  hasS <- any(mn$set == "S"); n_inv <- iv_main[["n"]] + iv_main[["na"]] + iv_cp[["n"]] + iv_cp[["na"]]
  EL[[length(EL) + 1]] <- data.frame(op = op, cand = cd, n_rows = nrow(mn), n_cases = length(unique(mn$key)), has_S = hasS, n_silent = n_sil, n_exc = n_exc, n_cov_fail = n_cov, n_fail_main = n_fail,
    inv_cells_main = iv_main[["n"]], inv_na_main = iv_main[["na"]], inv_cells_cp1 = iv_cp[["n"]], inv_na_cp1 = iv_cp[["na"]], n_cp1_rows = nrow(cp), n_cp1_fail = n_cpf, mut_ok = mut_ok,
    eligible = n_fail == 0 && n_sil == 0 && n_exc == 0 && n_cov == 0 && n_inv == 0 && n_cpf == 0 && mut_ok && hasS,
    eligible_noS = n_fail == 0 && n_sil == 0 && n_exc == 0 && n_cov == 0 && n_inv == 0 && n_cpf == 0 && mut_ok,
    guaranteed = cd %in% GUAR[[op]] && C6, t_S = tS(cd, op), rss_S = mS(cd, op), stringsAsFactors = FALSE)
}
EL <- do.call(rbind, EL); cat("\n== ELEGGIBILITÀ (ricalcolata)\n"); print(EL, row.names = FALSE)
write.csv(EL, file.path(REV, "recalc_eligibility.csv"), row.names = FALSE)
# ---- regola di decisione (punti 2-5) -----------------------------------------------------------------------------------
dec <- function(op) { e <- EL[EL$op == op & EL$eligible, ]; gq <- e[e$guaranteed, ]; cl <- if (nrow(gq)) "garantiti" else "non garantiti"; d <- if (nrow(gq)) gq else e
  if (!nrow(d)) return(data.frame(op = op, class = "nessuno", members = "", chosen = NA, reason = "nessun eleggibile"))
  if ("F1" %in% d$cand) { f1 <- d[d$cand == "F1", ]; o <- d[d$cand != "F1", ]
    ok <- !is.na(o$t_S) & !is.na(o$rss_S) & o$t_S <= f1$t_S / 2 & o$rss_S <= 1.5 * f1$rss_S
    if (any(ok)) { w <- o[ok, ][which.min(o$t_S[ok]), ]; return(data.frame(op = op, class = cl, members = paste(d$cand, collapse = ","), chosen = w$cand, reason = "sostituisce F1")) }
    return(data.frame(op = op, class = cl, members = paste(d$cand, collapse = ","), chosen = "F1", reason = sprintf("F1 resta (altri della classe: %s; tempi %s s vs F1 %.1f s)",
      paste(o$cand, collapse = ","), paste(sprintf("%.1f", o$t_S), collapse = ","), f1$t_S))) }
  w <- d[which.min(d$t_S), ]; data.frame(op = op, class = cl, members = paste(d$cand, collapse = ","), chosen = w$cand, reason = "il più veloce") }
DEC <- do.call(rbind, lapply(c("O1", "O2", "O3"), dec)); cat("\n== DECISIONE (ricalcolata)\n"); print(DEC, row.names = FALSE)
cat("un solo helper: scelte", paste(DEC$chosen, collapse = ","), "-> unico" , length(unique(DEC$chosen)) == 1, "\n")
# sensibilità: F5 in O1 senza il requisito S
cat("sensibilità: eleggibili O1 senza requisito S:", paste(EL$cand[EL$op == "O1" & EL$eligible_noS], collapse = ","), "\n")
write.csv(DEC, file.path(REV, "recalc_decision.csv"), row.names = FALSE)

# ---- CP-3 su H (F1) ---------------------------------------------------------------------------------------------------
h1 <- O1[O1$set == "H" & O1$cand == "F1", ]; h2 <- O2[O2$set == "H" & O2$cand == "F1" & O2$op %in% c("O2s2", "O2s8"), ]; h3 <- O2[O2$set == "H" & O2$cand == "F1" & O2$op == "O3", ]
cat("\n== CP-3 (F1 su H): O1", sum(h1$fail_eff), "/", nrow(h1), " (silenziosi", sum(h1$silent), ", eccezioni", sum(h1$exc > 0), ")  O2", sum(h2$fail), "/", nrow(h2), "  O3", sum(h3$fail), "/", nrow(h3),
    " | max ratio O2", format(max(h2$max_ratio), digits = 3), " O3", format(max(h3$max_ratio), digits = 3), " | casi H unici O1", length(unique(h1$key)), "\n")
cat("H: altri candidati (descrittivo) O1 fail:", paste(sapply(CANDS$O1, function(cd) { d <- O1[O1$set == "H" & O1$cand == cd & O1$applicable %in% TRUE, ]; sprintf("%s %d/%d", cd, sum(d$fail_eff), nrow(d)) }), collapse = "; "), "\n")
cat("H: copertura V0 invalida", sum(CV$v0[CV$set == "H"] %in% FALSE), "/", sum(CV$set == "H"), "; F4 valida", sum(CV$f4[CV$set == "H"] %in% TRUE), "\n")

# ---- numeri chiave -----------------------------------------------------------------------------------------------------
v0 <- O1[O1$set == "N" & O1$cand == "V0", ]
cat("\n== V0 in N: casi", nrow(v0), "; errore silenzioso", sum(v0$silent), "; unione", sum(v0$uv), "; diff. sim.", sum(v0$sv), "; con coppie", sum(v0$pv_any > 0), "; con celle", sum(v0$cell_viol > 0),
    "; unione E coppie", sum(v0$uv & v0$pv_any > 0), "; solo unione", sum(v0$uv & !(v0$pv_any > 0)), "; somma coppie violate (ordinate any)", sum(v0$pv_any), "; eccezioni", sum(v0$exc > 0), "\n")
cvN <- CV[CV$set == "N", ]; cat("copertura N: V0 invalida", sum(cvN$v0 %in% FALSE), "/", nrow(cvN), sprintf("(%.1f %%)", 100 * mean(cvN$v0 %in% FALSE)), "; F4 valida", sum(cvN$f4 %in% TRUE), "/", nrow(cvN), "\n")
f2N <- O1[O1$set == "N" & O1$cand == "F2", ]; cat("F2 in N: celle invalide dopo l'arrotondamento", sum(f2N$n_inv_rounded, na.rm = TRUE), "in", sum(f2N$n_inv_rounded > 0, na.rm = TRUE), "nulli; ancora invalide", sum(f2N$still_invalid, na.rm = TRUE), "; NA", sum(is.na(f2N$still_invalid)), "; errori silenziosi F2 su quei nulli", sum(f2N$silent[f2N$n_inv_rounded > 0]), "\n")
print(f2N[f2N$n_inv_rounded > 0, c("key", "n_inv_rounded", "still_invalid", "silent", "cell_viol")])
cat("F1/F3 celle invalide (validity_F1F3.csv) per insieme:\n"); VF$set <- subset_of(VF$key); print(aggregate(cbind(n_invalid_F1F3, n_empty_F1F3) ~ set, VF, sum))
cat("F4 celle invalide (c5) per insieme:\n"); print(aggregate(n_invalid ~ set, C5[!is.na(C5$n_invalid), ], sum)); print(C5[C5$n_invalid %in% 1:1e9, c("key", "n_invalid")])
# C-0
ns <- read.csv("results/S1.3/S1.3_null_summary.csv"); ns$key <- sprintf("%s_%s_%s_%02d", ns$archetype, ns$roi_id, ns$model, ns$rep)
c0 <- do.call(rbind, lapply(KN, function(k) { m <- meta[[k]]; data.frame(key = k, n = m$n, c1 = m$c0$c1, ov = m$c0$c2_overlap, sd = m$c0$c2_symdiff) }))
c0 <- merge(c0, ns[, c("key", "n", "c1", "c2_overlap", "c2_symdiff")], by = "key", suffixes = c("", ".r"), all.x = TRUE)
rl <- function(x, y) ifelse(x == y, 0, abs(x - y) / pmax(abs(y), 1e-300))
c0$ok <- c0$n == c0$n.r & rl(c0$c1, c0$c1.r) <= 1e-12 & rl(c0$ov, c0$c2_overlap) <= 1e-12 & rl(c0$sd, c0$c2_symdiff) <= 1e-12
cat("\nC-0:", sum(c0$ok %in% TRUE), "/", nrow(c0), "(righe S1.3 senza corrispondenza:", sum(is.na(c0$n.r)), "; n uguale:", sum(c0$n == c0$n.r, na.rm = TRUE), "; max rel c1/ov/sd:",
    format(max(rl(c0$c1, c0$c1.r)), digits = 3), format(max(rl(c0$ov, c0$c2_overlap)), digits = 3), format(max(rl(c0$sd, c0$c2_symdiff)), digits = 3), "; identici bit per bit c1/ov/sd:", sum(c0$c1 == c0$c1.r), sum(c0$ov == c0$c2_overlap), sum(c0$sd == c0$c2_symdiff), ")\n")
# C-5
C5$ok <- C5$n_over %in% 0 & C5$same_lost %in% TRUE & C5$same_gained %in% TRUE
cat("C-5 per insieme (ok / casi; parità; ambigui; errori):\n"); print(do.call(rbind, lapply(split(C5, C5$set), function(d) data.frame(set = d$set[1], ok = sum(d$ok), n = nrow(d), ties = sum(d$n_tie, na.rm = TRUE), cases_tie = sum(d$n_tie > 0, na.rm = TRUE), amb = sum(d$n_ambiguous, na.rm = TRUE), err = sum(!is.na(d$error))))))
nh <- C5$set != "H"; cat("C-5 senza H:", sum(C5$ok[nh]), "/", sum(nh), "; con H:", sum(C5$ok), "/", nrow(C5), "; casi con generatori attesi (meta has_pts) senza H:", sum(M$has_pts & M$set != "H"), "\n")
# CP-1
cat("\nCP-1 (O1, FAIL su input validi / validi; esclusi):\n")
for (cd in CANDS$O1) { d <- O1[O1$cand == cd & O1$set %in% CPS & O1$applicable %in% TRUE, ]; cat(sprintf("  %s: %d/%d; esclusi %d (T %d/%d, P %d/%d)\n", cd, sum(d$fail_eff[d$valid_input %in% TRUE]), sum(d$valid_input %in% TRUE), sum(!d$valid_input %in% TRUE),
  sum(d$fail_eff[d$valid_input %in% TRUE & d$set == "CP1T"]), sum(d$valid_input %in% TRUE & d$set == "CP1T"), sum(d$fail_eff[d$valid_input %in% TRUE & d$set == "CP1P"]), sum(d$valid_input %in% TRUE & d$set == "CP1P"))) }
cat("CP-1 (O2/O3, casi con FAIL su input validi):\n")
for (cd in CANDS$O2) for (op in c("O2s2", "O2s8", "O3")) { d <- O2[O2$cand == cd & O2$set %in% CPS & O2$op == op, ]; if (nrow(d)) cat(sprintf("  %s %s: %d/%d\n", cd, op, sum(d$fail[d$valid_input %in% TRUE]), sum(d$valid_input %in% TRUE))) }
cat("CP-1 input invalidi per sottoinsieme:", paste(names(table(M$set[M$set %in% CPS & !M$valid_input])), table(M$set[M$set %in% CPS & !M$valid_input]), collapse = " "), "\n")
# B-1
ref <- read.csv("results/S1.3/S1.3_verdicts.csv", stringsAsFactors = FALSE)
cat("\nB-1: righe S1.3_verdicts", nrow(ref), "; sezione B", sum(ref$section == "B"), "; B con confirmed non NA", sum(ref$section == "B" & !is.na(ref$confirmed)), "\n")
for (cd in c("V0", "F1", "F2", "F3", "F4", "F5")) { x <- read.csv(file.path(D, "b1", sprintf("verdicts_%s.csv", cd)), stringsAsFactors = FALSE)
  same_rows <- nrow(x) == nrow(ref) && all(x$section == ref$section & x$id == ref$id)
  kb <- ref$section == "B"; eqv <- x$value == ref$value; eqo <- x$outcome == ref$outcome; eqc <- (x$confirmed == ref$confirmed) | (is.na(x$confirmed) & is.na(ref$confirmed))
  cat(sprintf("  %s: stesse righe %s; B value uguali %d/%d; B outcome uguali %d/%d; B confirmed uguali %d/%d; righe non-B diverse (value) %d\n", cd, same_rows, sum(eqv[kb]), sum(kb), sum(eqo[kb]), sum(kb), sum(eqc[kb], na.rm = TRUE), sum(kb), sum(!eqv[!kb])))
  if (any(!eqv[kb])) print(cbind(ref[kb & !eqv, c("id", "value")], cand = x$value[kb & !eqv])) }
# K: i 26 casi noti
cat("\n== K\n"); ns$fail <- ns$c2_overlap > 1e-9 | ns$c2_symdiff > 1e-9; uf <- ns$key[ns$fail]; cat("unioni fallite in S1.3 (c2 > 1e-9):", length(uf), "; in K", sum(uf %in% KK), "; V0 union_viol in quei nulli:", sum(v0$uv[match(uf, v0$key)]), "\n")
p1 <- read.csv("results/S1.3b/S1.3b_posthoc_P1_pairs_v0.csv"); p1$key <- sprintf("%s_%s_%s_%02d", p1$archetype, p1$roi_id, p1$model, p1$rep)
pk <- rbind(p1[, c("key", "i", "j")], data.frame(key = "A5_r5_RSA_15", i = 71, j = 175))
pk$inK <- pk$key %in% KK
pk$v0_viol <- vapply(seq_len(nrow(pk)), function(q) { k <- pk$key[q]; o <- readRDS(file.path(D, "res", paste0(k, if (file.exists(file.path(D, "res", paste0(k, "_o12.rds")))) "_o12" else "_o1", ".rds")))
  v <- o$viol$o1_V0; if (is.null(v)) return(FALSE); a <- min(pk$i[q], pk$j[q]); b <- max(pk$i[q], pk$j[q]); any(v$what %in% c("pair_ij", "pair_ji") & pmin(v$i, v$j) == a & pmax(v$i, v$j) == b) }, TRUE)
print(pk); cat("casi noti:", length(uf) + nrow(pk), "; nulli distinti", length(unique(c(uf, pk$key))), "; K == unione:", setequal(KK, unique(c(uf, pk$key))), "\n")
k_v0 <- O1[O1$key %in% KK & O1$cand == "V0", ]; cat("V0 silent in K:", sum(k_v0$silent), "/", nrow(k_v0), "\n")
# ---- mtime: il riconfronto D-3 (compare.R modificato 2026-10-15 02:04:04) ha toccato solo S? -------------------------
tfix <- file.mtime("tools/S1.3d_compare.R"); MT$set <- subset_of(MT$key)
cat("\n== mtime res rispetto a compare.R (", format(tfix), "): file dopo la modifica per insieme:\n"); print(table(MT$set, MT$mtime > tfix))
sfi <- file.info(file.path(D, "res", sf_)); print(data.frame(file = sf_, mtime = format(sfi$mtime), after_fix = sfi$mtime > tfix), row.names = FALSE)
cat("\nFINE\n")
