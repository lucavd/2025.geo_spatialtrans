.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

OUTD <- "results/S1.4b/review/recalc"
FF <- read.csv("results/S1.4b/S1.4b_folds.csv", stringsAsFactors = FALSE)
F <- read.csv(file.path(OUTD, "recalc_folds.csv"), stringsAsFactors = FALSE)
F <- F[F$variant != "V-b", ]; F$variant[F$variant == "V-a"] <- "V-a/V-b"
W <- reshape(F[, c("archetype","fold","variant","metric","dist","verdict")], idvar = c("archetype","fold","variant"), timevar = "metric", direction = "wide")
ex <- unique(F[, c("archetype","fold","variant","trunc_all_med","loss_all_med","b7_sim_med","b7_real","yard_med","yard_max","metric")])
ex <- ex[ex$metric == "b7", ]
M <- merge(merge(FF, W, by = c("archetype","fold","variant")), ex, by = c("archetype","fold","variant"))
cat("righe:", nrow(FF), nrow(M), "\n")
d <- c(ks_nc = max(abs(M$ks_nc - M$dist.ks_nc)), ks_an = max(abs(M$ks_an - M$dist.ks_an)), d_sp = max(abs(M$d_sp - M$dist.d_sp)),
       b7 = max(abs(M$b7 - M$dist.b7)), trunc = max(abs(M$trunc_all - M$trunc_all_med)), loss = max(abs(M$loss_all - M$loss_all_med)),
       b7_sim = max(abs(M$b7_sim - M$b7_sim_med)), b7_real = max(abs(M$b7_real.x - M$b7_real.y)),
       y_med = max(abs(M$y_b7_med - M$yard_med)), y_max = max(abs(M$y_b7_max - M$yard_max)))
print(signif(d, 3))
cat("verdetti fold identici:", sum(M$v_b1 == M$verdict.ks_nc & M$v_b2 == M$verdict.ks_an & M$v_b3 == M$verdict.d_sp & M$v_b7 == M$verdict.b7), "/", nrow(M), "\n")
# margini: fold piu' vicini a una soglia (rischio di ribaltamento)
F2 <- read.csv(file.path(OUTD, "recalc_folds.csv"), stringsAsFactors = FALSE); F2 <- F2[F2$variant != "V-b", ]
F2$rel_to_med <- (F2$dist - F2$yard_med) / F2$yard_med; F2$rel_to_max <- (F2$dist - F2$yard_max) / F2$yard_max
F2$near <- pmin(abs(F2$rel_to_med), abs(F2$rel_to_max))
o <- F2[order(F2$near), c("archetype","fold","variant","metric","dist","yard_med","yard_max","verdict","near")]
print(head(o, 12), row.names = FALSE, digits = 4)
