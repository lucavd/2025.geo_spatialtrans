.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
TB <- readRDS("/mnt/micron/geo_spatialtrans/S1.3d/review_recalc/recalc_tables.rds"); O1 <- TB$O1; O2 <- TB$O2
d <- O2[O2$cand == "F5" & O2$fail, ]; print(table(d$set, d$op)); print(d[d$set %in% c("N","U1","U3","U4","U5","S"), c("key","op","n_viol","n_missing","n_spurious","max_ratio")])
cat("V0 O1 silent per insieme:\n"); v <- O1[O1$cand == "V0" & O1$applicable %in% TRUE, ]; print(table(v$set, v$silent))
cat("F4 O1 applicabile per insieme:\n"); print(table(O1$set[O1$cand=="F4" & O1$applicable %in% TRUE]))
cat("max ratio O1 F1 coppie / celle per insieme:\n"); f <- O1[O1$cand == "F1", ]; print(aggregate(cbind(pmax_ratio, cmax_ratio) ~ set, f, max))
cat("max ratio O2 F1 per insieme e op:\n"); f <- O2[O2$cand == "F1", ]; print(aggregate(max_ratio ~ set + op, f, max))
cat("U3 F1 (cella vuota dopo set_precision):\n"); print(O1[O1$set=="U3" & O1$cand=="F1", c("key","silent","cell_viol","cmax_ratio")])
