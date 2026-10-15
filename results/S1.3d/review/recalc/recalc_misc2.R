.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
TB <- readRDS("/mnt/micron/geo_spatialtrans/S1.3d/review_recalc/recalc_tables.rds"); O1 <- TB$O1; CV <- TB$CV
v <- O1[O1$cand == "V0" & O1$set == "N", ]; v$cov <- CV$v0[match(v$key, CV$key)]
cat("V0 N: silenziosi", sum(v$silent), "; di cui copertura valida", sum(v$silent & v$cov %in% TRUE), "; copertura invalida senza errore", sum(!v$silent & v$cov %in% FALSE), "\n")
