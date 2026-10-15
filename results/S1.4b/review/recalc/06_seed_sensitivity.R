.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

# recalc — sensibilita' della regola di decisione alla scelta dei 10 semi (bootstrap sui semi e leave-one-seed-out)
OUTD <- "results/S1.4b/review/recalc"
fs <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.4b/loo", full.names = TRUE))
L <- lapply(fs, readRDS)
vd <- function(x, med, mx) ifelse(x <= med, 2L, ifelse(x <= mx, 1L, 0L))  # 2 PASS 1 WARN 0 FAIL
arch_rule <- function(v) { if (sum(v == 2L) >= 4) 2L else if (sum(v == 0L) >= 2) 0L else 1L }
mets <- c("ks_nc", "ks_an", "d_sp", "b7")
# array per fold: valori per seme [10, metric, variant]
prep <- lapply(L, function(x) { s <- x$m[x$m$seed_idx %in% 1:10, ]
  arr <- array(NA, c(10, 4, 2), dimnames = list(NULL, mets, c("V-c", "V-a")))
  for (v in c("V-c", "V-a")) { sv <- s[s$variant == v, ]; sv <- sv[order(sv$seed_idx), ]
    arr[, , v] <- cbind(sv$ks_nc, sv$ks_an, sv$d_sp, abs(sv$b7 - x$b7_real)) }
  list(A = x$archetype, k = x$fold, arr = arr,
       med = c(apply(x$yard[, 2:4], 2, median), median(x$yard_b7)), mx = c(apply(x$yard[, 2:4], 2, max), max(x$yard_b7))) })
decide <- function(idx_list) {   # idx_list: per fold, indici dei semi
  cnt <- c("V-c" = 0, "V-a" = 0); per <- list()
  for (A in paste0("A", 1:6)) { ii <- which(vapply(prep, function(p) p$A, "") == A)
    for (v in c("V-c", "V-a")) for (mi in 1:4) {
      fv <- vapply(ii, function(f) { p <- prep[[f]]; vd(median(p$arr[idx_list[[f]], mi, v]), p$med[mi], p$mx[mi]) }, 1L)
      o <- arch_rule(fv); per[[paste(A, v, mets[mi])]] <- o; if (o == 2L) cnt[v] <- cnt[v] + 1 } }
  list(cnt = cnt, per = unlist(per)) }
base <- decide(rep(list(1:10), length(prep))); cat("base:", base$cnt, "\n")
set.seed(20261015); B <- 2000
bt <- t(vapply(1:B, function(b) { idx <- lapply(seq_along(prep), function(f) sample.int(10, 10, TRUE)); d <- decide(idx); c(d$cnt, d$per) }, numeric(2 + length(base$per))))
cnt_ab <- bt[, 2]; cnt_c <- bt[, 1]
cat(sprintf("bootstrap semi (B=%d, ricampionamento indipendente per fold): P(ab>c)=%.3f P(tie)=%.3f P(c>ab)=%.3f\n", B, mean(cnt_ab > cnt_c), mean(cnt_ab == cnt_c), mean(cnt_ab < cnt_c)))
print(table(ab = cnt_ab, c = cnt_c))
stab <- colMeans(bt[, -(1:2)] == matrix(base$per, B, length(base$per), byrow = TRUE))
uns <- sort(stab[stab < 0.95]); print(round(uns, 3))
# leave-one-seed-out (stesso seme tolto in tutti i fold)
lo <- t(vapply(1:10, function(s) decide(rep(list(setdiff(1:10, s)), length(prep)))$cnt, numeric(2)))
print(cbind(seed_out = 1:10, lo))
# semi del fold A6 fold 3 V-c ks_an
p <- prep[[which(vapply(prep, function(p) p$A == "A6" && p$k == 3, TRUE))]]
cat("A6 fold3 V-c ks_an per seme:", round(p$arr[, "ks_an", "V-c"], 4), " soglia med", round(p$med[2], 5), " mediana", round(median(p$arr[, "ks_an", "V-c"]), 5), "\n")
write.csv(data.frame(cell = names(stab), p_same_as_base = stab), file.path(OUTD, "recalc_seed_bootstrap_stability.csv"), row.names = FALSE)
write.csv(data.frame(cnt_c = cnt_c, cnt_ab = cnt_ab), file.path(OUTD, "recalc_seed_bootstrap_counts.csv"), row.names = FALSE)
