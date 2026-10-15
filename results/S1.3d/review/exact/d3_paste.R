# review/exact/d3_paste.R — punto (e): quali interi double danno paste() diverso da sprintf("%d")? (R 4.6, scipen di default)
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(jsonlite) })
cat("scipen =", getOption("scipen"), "\n")
k <- as.double(0:2000000); bad <- k[paste(k) != sprintf("%d", as.integer(k))]
cat("interi in [0, 2e6] con paste(double) != sprintf: ", length(bad), "\n"); print(head(bad, 30))
cat("minimo:", min(bad), " tutti multipli di 1e5:", all(bad %% 1e5 == 0), "\n")
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; mf <- list.files(file.path(D, "geom"), pattern = "_meta\\.json$", full.names = TRUE)
M <- do.call(rbind, lapply(mf, function(f) { m <- fromJSON(f); data.frame(key = m$key, n = m$n, nx2 = m$gs$nx8 * 4, ny2 = m$gs$ny8 * 4) }))
M$S <- grepl("^S_", M$key)
cat("casi:", nrow(M), "; fuori da S: max n =", max(M$n[!M$S]), " max nx2 =", max(M$nx2[!M$S]), " max ny2 =", max(M$ny2[!M$S]), "\n")
cat("S: n =", M$n[M$S], " nx2 =", M$nx2[M$S], " ny2 =", M$ny2[M$S], "\n")
cat("fuori da S, casi con n >= 1e5 o nx2/ny2 >= 1e5:", sum(!M$S & (M$n >= 1e5 | M$nx2 >= 1e5 | M$ny2 >= 1e5)), "\n")
# chiavi di bin a 8 µm: indici < nx8; e indici ix/iy usati nei confronti a 2 µm < nx2
for (s in c(2, 8)) { T <- read_parquet(file.path(D, "truth", sprintf("S_cperf_o2_s%d.parquet", s)))
  sel <- T$i %in% c(1e5, 2e5, 3e5, 4e5); cat(sprintf("S, s = %d: pezzi esatti delle celle 1e5..4e5 = %d (-> %d violazioni fittizie se mancanti + spuri)\n", s, sum(sel), 2 * sum(sel)))
  cat("   per cella:", paste(names(table(T$i[sel])), table(T$i[sel]), collapse = "; "), "\n")
  kd <- paste(as.double(T$i), as.double(T$ix), as.double(T$iy)); ki <- sprintf("%d %d %d", T$i, T$ix, T$iy)
  cat("   pezzi di S con chiave paste(double) diversa:", sum(kd != ki), " di cui fuori dalle celle 1e5..4e5:", sum(kd != ki & !sel), "\n") }
