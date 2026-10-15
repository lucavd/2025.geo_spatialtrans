# results/S1.3d/review/recalc/recalc_validity.R — revisione «ricalcolo» S1.3d: validità degli ingressi (st_is_valid sf, GEOS 3.12) e
# delle celle F1/F3 (geos_set_precision g = 1e-6, libgeos del pacchetto geos) per TUTTE le chiavi (N, H, U, CP1, S), con lettore parquet proprio;
# controllo dei perimetri P di meta (ricalcolati). Scrive /mnt/micron/geo_spatialtrans/S1.3d/review_recalc/input_validity.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow); library(jsonlite); library(parallel) })
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; HEAVY <- file.path(D, "review_recalc")
read_cells <- function(f, n) { t <- as.data.frame(read_parquet(f)); t <- t[order(t$cell, t$part, t$ring, t$k), ]
  by <- split(t, t$cell); out <- vector("list", n)
  for (q in seq_len(n)) { d <- by[[as.character(q)]]
    if (is.null(d)) { out[[q]] <- st_polygon(); next }
    parts <- lapply(split(d, d$part), function(dp) lapply(split(dp, dp$ring), function(dr) { m <- cbind(dr$x, dr$y); rbind(m, m[1, ]) }))
    out[[q]] <- if (length(parts) == 1) st_polygon(parts[[1]]) else st_multipolygon(unname(parts)) }
  st_sfc(out) }
keys <- c(readLines(file.path(D, "keys_N.txt")), readLines(file.path(D, "keys_H.txt")), readLines(file.path(D, "keys_U.txt")), readLines(file.path(D, "keys_CP1.txt")), "S_cperf")
R <- mclapply(keys, function(k) { m <- fromJSON(file.path(D, "geom", paste0(k, "_meta.json")))
  t <- read_cells(file.path(D, "geom", paste0(k, "_cells.parquet")), m$n)
  v <- st_is_valid(t); P <- as.numeric(st_length(st_cast(st_boundary(t), "MULTILINESTRING")))
  s <- geos_set_precision(as_geos_geometry(t), 1e-6); vs <- geos_is_valid(s)
  data.frame(key = k, n = m$n, n_invalid_input = sum(!v | is.na(v)), n_empty_input = sum(st_is_empty(t)), n_invalid_F1F3 = sum(!vs | is.na(vs)), n_empty_F1F3 = sum(geos_is_empty(s)),
             P_maxrel = max(abs(P - m$P) / m$P), stringsAsFactors = FALSE) }, mc.cores = 28, mc.preschedule = FALSE)
err <- vapply(R, function(r) !is.data.frame(r), TRUE); if (any(err)) { print(R[err][1:3]); stop("errori") }
R <- do.call(rbind, R); write.csv(R, file.path(HEAVY, "input_validity.csv"), row.names = FALSE)
st <- ifelse(R$key %in% readLines(file.path(D, "keys_N.txt")), "N", ifelse(R$key %in% readLines(file.path(D, "keys_H.txt")), "H", ifelse(R$key == "S_cperf", "S", substr(R$key, 1, ifelse(startsWith(R$key, "CP1"), 4, 2)))))
print(aggregate(cbind(n_invalid_input, n_empty_input, n_invalid_F1F3, n_empty_F1F3) ~ st, R, sum)); print(aggregate(cbind(cases_inv_input = n_invalid_input > 0) ~ st, R, sum))
cat("perimetri: max rel |P_ricalcolato - P_meta| / P =", format(max(R$P_maxrel), digits = 3), "\n")
VF <- read.csv(file.path(D, "checks", "validity_F1F3.csv")); m <- merge(VF, R, by = "key", suffixes = c(".off", ".rec"))
cat("validity_F1F3.csv vs ricalcolo: chiavi in comune", nrow(m), "; n_invalid diverso", sum(m$n_invalid_F1F3.off != m$n_invalid_F1F3.rec), "; n_empty diverso", sum(m$n_empty_F1F3.off != m$n_empty_F1F3.rec), "\n")
cat("FINE\n")
