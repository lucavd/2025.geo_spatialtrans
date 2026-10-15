.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

# recalc — ricalcolo indipendente dei verdetti B, CP-1, previsioni e regola (S1.4b), dai grezzi loo/*.rds e cp1/*.csv
D <- "/mnt/micron/geo_spatialtrans/S1.4b"; OUTD <- "results/S1.4b/review/recalc"
fs <- sort(list.files(file.path(D, "loo"), pattern = "\\.rds$", full.names = TRUE))
stopifnot(length(fs) == 30)
vd <- function(x, med, mx) ifelse(x <= med, "PASS", ifelse(x <= mx, "WARN", "FAIL"))
fold_rows <- list(); chk <- list()
for (f in fs) {
  x <- readRDS(f); m <- x$m
  stopifnot(x$mutant == "none")
  s <- m[m$seed_idx %in% 1:10, ]
  # V-a == V-b sulle metriche B?
  a <- s[s$variant == "V-a", ]; b <- s[s$variant == "V-b", ]
  a <- a[order(a$seed_idx), ]; b <- b[order(b$seed_idx), ]
  chk[[f]] <- data.frame(file = basename(f), ab_identical = isTRUE(all.equal(a[, -2], b[, -2], check.attributes = FALSE, tolerance = 0)),
                         ab_maxdiff = max(abs(as.matrix(a[, 3:9]) - as.matrix(b[, 3:9]))), nseeds = length(unique(s$seed_idx)),
                         yard_n = nrow(x$yard), yardb7_n = length(x$yard_b7))
  for (v in c("V-c", "V-a", "V-b")) {
    sv <- s[s$variant == v, ]; stopifnot(nrow(sv) == 10, setequal(sv$seed_idx, 1:10))
    d <- c(ks_nc = median(sv$ks_nc), ks_an = median(sv$ks_an), d_sp = median(sv$d_sp),
           b7 = median(abs(sv$b7 - x$b7_real)))
    d_b7_alt <- abs(median(sv$b7) - x$b7_real)
    med <- c(apply(x$yard[, c("ks_nc", "ks_an", "d_sp")], 2, median), b7 = median(x$yard_b7))
    mx  <- c(apply(x$yard[, c("ks_nc", "ks_an", "d_sp")], 2, max), b7 = max(x$yard_b7))
    fold_rows[[length(fold_rows) + 1]] <- data.frame(archetype = x$archetype, method = x$method, fold = x$fold, roi_test = x$roi_test,
      variant = v, metric = names(d), dist = unname(d), yard_med = unname(med), yard_max = unname(mx),
      verdict = vd(unname(d), unname(med), unname(mx)),
      verdict_b7_alt = ifelse(names(d) == "b7", vd(d_b7_alt, med[["b7"]], mx[["b7"]]), NA),
      b7_sim_med = median(sv$b7), b7_real = x$b7_real,
      trunc_all_med = median(sv$trunc_all), loss_all_med = median(sv$loss_all),
      b8_circle = x$b8[["circle"]], b8_rule = x$b8[["rule"]],
      dat_med = median(x$checks$dat_med[x$checks$seed_idx %in% 1:10]), dat_real_med = median(x$checks$dat_real_med[x$checks$seed_idx %in% 1:10]),
      cv_vor = median(x$checks$cv_vor[x$checks$seed_idx %in% 1:10]), cv_lag = median(x$checks$cv_lag[x$checks$seed_idx %in% 1:10]),
      cv_lag_real = median(x$checks$cv_lag_real[x$checks$seed_idx %in% 1:10]),
      trunc_pair = median(x$checks$trunc_pair[x$checks$seed_idx %in% 1:10]), trunc_reg = median(x$checks$trunc_reg[x$checks$seed_idx %in% 1:10]),
      stringsAsFactors = FALSE)
  }
}
F <- do.call(rbind, fold_rows); CK <- do.call(rbind, chk)
write.csv(F, file.path(OUTD, "recalc_folds.csv"), row.names = FALSE)
write.csv(CK, file.path(OUTD, "recalc_ab_identity.csv"), row.names = FALSE)
cat("V-a == V-b identici in tutti i fold:", all(CK$ab_identical), " max diff:", max(CK$ab_maxdiff), "\n")
cat("yard n:", unique(CK$yard_n), unique(CK$yardb7_n), "\n")
arch_rule <- function(v) { np <- sum(v == "PASS"); nf <- sum(v == "FAIL"); if (np >= 4) "PASS" else if (nf >= 2) "FAIL" else "WARN" }
mname <- c(ks_nc = "B-1", ks_an = "B-2", d_sp = "B-3", b7 = "B-7")
AR <- do.call(rbind, lapply(split(F, list(F$archetype, F$variant, F$metric), drop = TRUE), function(g) {
  g <- g[order(g$fold), ]
  data.frame(archetype = g$archetype[1], variant = g$variant[1], id = mname[[g$metric[1]]],
             folds = paste(substr(g$verdict, 1, 1), collapse = ""), dists = paste(signif(g$dist, 2), collapse = " "),
             outcome = arch_rule(g$verdict),
             outcome_b7_alt = if (g$metric[1] == "b7") arch_rule(g$verdict_b7_alt) else NA,
             b7_sim_mean = mean(g$b7_sim_med), b7_real_mean = mean(g$b7_real), stringsAsFactors = FALSE) }))
AR <- AR[order(AR$variant, AR$id, AR$archetype), ]
# B-4
B4 <- do.call(rbind, lapply(split(F[F$metric == "ks_nc", ], list(F$archetype[F$metric == "ks_nc"], F$variant[F$metric == "ks_nc"]), drop = TRUE), function(g)
  data.frame(archetype = g$archetype[1], variant = g$variant[1], b4 = mean(g$trunc_all_med), loss = mean(g$loss_all_med),
             outcome = ifelse(mean(g$trunc_all_med) <= 0.05, "PASS", "FAIL"),
             b8_circle = mean(g$b8_circle), b8_rule = mean(g$b8_rule), b9_lag = median(g$dat_med), b9_lag_real = median(g$dat_real_med),
             b9_mean = mean(g$dat_med), cv_vor = median(g$cv_vor), cv_lag = median(g$cv_lag), cv_lag_real = median(g$cv_lag_real),
             trunc_pair = mean(g$trunc_pair), trunc_reg = mean(g$trunc_reg))))
B4 <- B4[order(B4$variant, B4$archetype), ]
write.csv(AR, file.path(OUTD, "recalc_archetype_B.csv"), row.names = FALSE)
write.csv(B4, file.path(OUTD, "recalc_B4_B8_B9.csv"), row.names = FALSE)
options(width = 200)
print(AR[AR$variant != "V-b", ], row.names = FALSE)
print(B4[B4$variant != "V-b", ], row.names = FALSE, digits = 4)
# ---- CP-1
cf <- sort(list.files(file.path(D, "cp1"), pattern = "\\.csv$", full.names = TRUE)); stopifnot(length(cf) == 30)
CP <- do.call(rbind, lapply(cf, function(f) { d <- read.csv(f)
  p <- d[grepl("^perm", d$cond), ]; stopifnot(nrow(p) == 20)
  data.frame(archetype = d$archetype[1], roi = d$roi[1], method = d$method[1], conds = paste(d$cond, collapse = ";"),
             m1_vor = d$m1[d$cond == "vor"], m1_lag = d$m1[d$cond == "lag"], m1_perm_min = min(p$m1), m1_perm_mean = mean(p$m1),
             m2_vor = d$m2[d$cond == "vor"], m2_lag = d$m2[d$cond == "lag"], m2_perm_min = min(p$m2), m2_perm_mean = mean(p$m2),
             beats = (d$m1[d$cond == "lag"] < min(p$m1)) & (d$m2[d$cond == "lag"] < min(p$m2)),
             beats_m1 = d$m1[d$cond == "lag"] < min(p$m1), beats_m2 = d$m2[d$cond == "lag"] < min(p$m2), stringsAsFactors = FALSE) }))
write.csv(CP, file.path(OUTD, "recalc_cp1_roi.csv"), row.names = FALSE)
cat("cond set:", unique(CP$conds), "\n")
CPA <- do.call(rbind, lapply(split(CP, CP$archetype), function(g) data.frame(archetype = g$archetype[1], n_beat = sum(g$beats),
  beats = sum(g$beats) >= 4, m1_vor = mean(g$m1_vor), m1_lag = mean(g$m1_lag), m1_perm = mean(g$m1_perm_mean),
  m2_vor = mean(g$m2_vor), m2_lag = mean(g$m2_lag), m2_perm = mean(g$m2_perm_mean),
  red_mean_ratio = mean(1 - g$m2_lag / g$m2_vor), red_ratio_means = 1 - mean(g$m2_lag) / mean(g$m2_vor),
  margin_m2_min = min(g$m2_perm_min - g$m2_lag), margin_m1_min = min(g$m1_perm_min - g$m1_lag))))
write.csv(CPA, file.path(OUTD, "recalc_cp1_archetype.csv"), row.names = FALSE)
print(CPA, row.names = FALSE, digits = 3)
n_arch_beat <- sum(CPA$beats); a_refuted <- n_arch_beat < 5
cat("CP-1 archetipi che battono:", n_arch_beat, "/6 ; (a) smentita:", a_refuted, "\n")
# ---- previsioni
g <- function(v, id, A) AR$outcome[AR$variant == v & AR$id == id & AR$archetype == A]
A14 <- paste0("A", 1:4); A15 <- paste0("A", 1:5); A16 <- paste0("A", 1:6)
P1 <- all(CPA$red_mean_ratio[CPA$archetype %in% A14] >= 0.5) && n_arch_beat == 6
b4 <- function(v, A) B4$b4[B4$variant == v & B4$archetype == A]
P2 <- all(sapply(A14, function(A) b4("V-a", A) <= 0.5 * b4("V-c", A)))
p3n <- sum(sapply(A16, function(A) g("V-a", "B-1", A) == "PASS" && g("V-a", "B-2", A) == "PASS"))
p3a <- sum(sapply(A16, function(A) g("V-a", "B-1", A) == "PASS")); p3b <- sum(sapply(A16, function(A) g("V-a", "B-2", A) == "PASS"))
P3 <- p3n >= 4
p4a <- sum(sapply(A16, function(A) g("V-a", "B-7", A) %in% c("PASS", "WARN"))); p4c <- sum(sapply(A15, function(A) g("V-c", "B-7", A) == "FAIL"))
P4 <- p4a >= 5 && p4c == 5
P5 <- all(B4$b9_lag[B4$variant == "V-a" & B4$archetype %in% A14] <= 0.10)
P6 <- all(sapply(c("V-c", "V-a", "V-b"), function(v) g(v, "B-1", "A5") != "PASS"))
cat(sprintf("P1 %s (red %s); P2 %s; P3 %s (%d; B1 %d | B2 %d); P4 %s (V-a/b %d/6; V-c FAIL %d/5); P5 %s; P6 %s\n",
  P1, paste(round(CPA$red_mean_ratio, 3), collapse = " "), P2, P3, p3n, p3a, p3b, P4, p4a, p4c, P5, P6))
cnt <- function(v) sum(AR$variant == v & AR$id %in% c("B-1", "B-2", "B-3", "B-7") & AR$outcome == "PASS")
cat("Conteggio PASS {B-1,B-2,B-3,B-7}: V-a", cnt("V-a"), " V-b", cnt("V-b"), " V-c", cnt("V-c"), "\n")
cnt_alt <- function(v) sum(AR$variant == v & ((AR$id %in% c("B-1", "B-2", "B-3") & AR$outcome == "PASS") | (AR$id == "B-7" & AR$outcome_b7_alt == "PASS")), na.rm = TRUE)
cat("Conteggio con B-7 alternativo |med(b7)-b7_real|: V-a", cnt_alt("V-a"), " V-c", cnt_alt("V-c"), "\n")
print(table(F$verdict[F$metric=="b7"], F$verdict_b7_alt[F$metric=="b7"]))
