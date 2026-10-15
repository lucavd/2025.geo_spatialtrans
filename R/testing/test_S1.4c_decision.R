# R/testing/test_S1.4c_decision.R — test sintetici delle funzioni di decisione di S1.4c (prima del run).
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S1.4c_decision.R
source("tools/S1.4c_decision.R")
res <- list(); ck <- function(name, ok) { res[[length(res) + 1]] <<- isTRUE(ok); cat(sprintf("%s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name)) }
ck("D1 fold_verdict PASS/WARN/FAIL", identical(fold_verdict(c(0.1, 0.3, 0.5), c(0.1, 0.2, 0.3, 0.4)), c("PASS", "WARN", "FAIL")))
ck("D1b fold_verdict al confine (= mediana -> PASS, = max -> WARN)", identical(fold_verdict(c(0.25, 0.4), c(0.1, 0.2, 0.3, 0.4)), c("PASS", "WARN")))
ck("D2 arch_verdict", identical(c(arch_verdict(c(rep("PASS", 4), "FAIL")), arch_verdict(c("PASS", "PASS", "PASS", "FAIL", "FAIL")),
                                  arch_verdict(c("PASS", "PASS", "PASS", "WARN", "FAIL"))), c("PASS", "FAIL", "WARN")))
# dati sintetici: 2 archetipi x 5 fold x 10 semi x 2 varianti x 1 metrica; K-x vicino (PASS), K-c lontano (FAIL)
yd <- list(); for (A in c("A1", "A2")) { yd[[A]] <- list(); for (k in 1:5) yd[[A]][[k]] <- list(B1 = c(0.1, 0.2, 0.3, 0.4)) }
d <- expand.grid(archetype = c("A1", "A2"), fold = 1:5, seed_idx = 1:10, variant = c("K-c", "K-1"), stringsAsFactors = FALSE)
d$B1 <- ifelse(d$variant == "K-1", 0.05, 0.9)
vt <- verdicts(d, yd, metrics = "B1")
ck("D3 verdicts: K-1 PASS e K-c FAIL in entrambi", all(vt$verdict[vt$variant == "K-1"] == "PASS") && all(vt$verdict[vt$variant == "K-c"] == "FAIL"))
sc <- score_of(vt, "B1"); ck("D4 score", sc[["K-1"]] == 2 && sc[["K-c"]] == 0)
# mediana sui semi: un seme anomalo non cambia il verdetto
d2 <- d; d2$B1[d2$variant == "K-1" & d2$seed_idx == 3] <- 9
ck("D5 mediana dei semi robusta a un seme anomalo", all(verdicts(d2, yd, "B1")$verdict[vt$variant == "K-1"] == "PASS"))
bs <- boot_scores(d, yd, "B1", B = 20)
ck("D6 bootstrap: forma e valori", all(dim(bs) == c(20, 2)) && all(bs[, "K-1"] == 2) && all(bs[, "K-c"] == 0))
dc <- decide(c("K-c" = 10, "K-1" = 13, "K-2" = 14), c("K-1" = 0.95, "K-2" = 0.85), c("K-1" = TRUE, "K-2" = TRUE))
ck("D7 decide: K-2 ha piu' punti ma p < 0.9 -> K-1", dc$choice == "K-1")
dc <- decide(c("K-c" = 10, "K-1" = 12, "K-2" = 14), c("K-1" = 0.99, "K-2" = 0.95), c("K-1" = TRUE, "K-2" = TRUE))
ck("D8 decide: K-1 margine 2 < 3 -> K-2", dc$choice == "K-2")
dc <- decide(c("K-c" = 10, "K-1" = 14, "K-2" = 14), c("K-1" = 0.99, "K-2" = 0.99), c("K-1" = TRUE, "K-2" = TRUE))
ck("D9 decide: parita' -> la piu' semplice (K-1)", dc$choice == "K-1")
dc <- decide(c("K-c" = 10, "K-1" = 14, "K-2" = 15), c("K-1" = 0.99, "K-2" = 0.99), c("K-1" = FALSE, "K-2" = FALSE))
ck("D10 decide: nessuna eleggibile -> K-c", dc$choice == "K-c")
ck("D11 auc: separazione perfetta = 1, inversa = 0, casuale ~ 0.5",
   auc(1:10, c(rep(FALSE, 5), rep(TRUE, 5))) == 1 && auc(10:1, c(rep(FALSE, 5), rep(TRUE, 5))) == 0 &&
   abs({ set.seed(1); auc(runif(4000), runif(4000) > 0.5) } - 0.5) < 0.03)
ck("D12 auc con pareggi = 0.5", auc(rep(1, 6), c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)) == 0.5)
yy <- rep(list(c(0.1, 0.2, 0.3, 0.4)), 5)
ck("D13 cp2_exceeds: 3/5 sopra il max -> TRUE; 2/5 -> FALSE", cp2_exceeds(c(0.5, 0.5, 0.5, 0.1, 0.1), yy) && !cp2_exceeds(c(0.5, 0.5, 0.4, 0.1, 0.1), yy))
nf <- sum(!unlist(res)); cat(sprintf("RESULT %d/%d PASS\n", length(res) - nf, length(res))); quit(status = if (nf) 1 else 0)
