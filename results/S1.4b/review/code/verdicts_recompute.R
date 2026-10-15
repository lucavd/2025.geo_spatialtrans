
# review/code/verdicts_recompute.R — revisore "code": verdetti B ricalcolati da zero (metriche di vc_poi.R) con
# implementazione propria della regola; V-c (banco), V-a/V-b, e V-c con l'operatore "Voronoi dei centri dei nuclei".
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
fs <- list.files("/tmp/code_vc", pattern = "fold[0-9]\\.rds$", full.names = TRUE)
cl <- function(v, y) ifelse(v <= median(y), "P", ifelse(v <= max(y), "W", "F"))
av <- function(o) if (sum(o == "P") >= 4) "PASS" else if (sum(o == "F") >= 2) "FAIL" else "WARN"
F <- do.call(rbind, lapply(fs, function(f) { o <- readRDS(f)
  do.call(rbind, lapply(unique(o$m$variant), function(v) { m <- o$m[o$m$variant == v, ]
    d <- c(ks_nc = median(m$ks_nc), ks_an = median(m$ks_an), d_sp = median(m$d_sp), b7 = median(abs(m$b7 - o$b7_real)))
    data.frame(A = o$A, k = o$k, variant = v, t(d), b7_sim = median(m$b7), b7_real = o$b7_real, trunc = median(m$trunc_all), loss = median(m$loss_all),
               v1 = cl(d[["ks_nc"]], o$yard$ks_nc), v2 = cl(d[["ks_an"]], o$yard$ks_an), v3 = cl(d[["d_sp"]], o$yard$d_sp), v7 = cl(d[["b7"]], o$yard_b7),
               stored_eq = if (all(is.na(m$stored_equal))) NA else all(m$stored_equal), at_p_vs_at = o$at_p_vs_at, poi_shift_med = o$poi_shift[["Median"]]) })) }))
F <- F[order(F$variant, F$A, F$k), ]
write.csv(F, "results/S1.4b/review/code/recomputed_folds.csv", row.names = FALSE)
Atab <- do.call(rbind, lapply(split(F, list(F$A, F$variant)), function(z) data.frame(A = z$A[1], variant = z$variant[1],
  B1 = av(z$v1), B2 = av(z$v2), B3 = av(z$v3), B7 = av(z$v7), f1 = paste(z$v1, collapse = ""), f2 = paste(z$v2, collapse = ""), f3 = paste(z$v3, collapse = ""), f7 = paste(z$v7, collapse = ""),
  B4 = mean(z$trunc), loss = mean(z$loss), b7_sim = mean(z$b7_sim), b7_real = mean(z$b7_real))))
Atab <- Atab[order(Atab$variant, Atab$A), ]; rownames(Atab) <- NULL
write.csv(Atab, "results/S1.4b/review/code/recomputed_archetypes.csv", row.names = FALSE)
print(Atab, digits = 4)
for (v in unique(Atab$variant)) cat(v, "score PASS {B1,B2,B3,B7} =", sum(unlist(Atab[Atab$variant == v, c("B1", "B2", "B3", "B7")]) == "PASS"), "\n")
cat("stored equality (tolerance 0) per fold:\n"); print(table(F$variant, F$stored_eq, useNA = "ifany"))
cat("POI-Voronoi vs generator-Voronoi, median |dA_T|/A_T per fold:\n"); print(tapply(F$at_p_vs_at[F$variant == "V-c_poiVor"], F$A[F$variant == "V-c_poiVor"], median))
