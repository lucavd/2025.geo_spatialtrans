# tools/S2.6_measures.R - S2.6: misure delle giunture (J-1, J-2, J-3) su un'esecuzione di tools/S2.6_pieces.R. Definizioni: docs/steps/S2.6.md.
# Blocchi di 8 um della griglia dell'area intera con il centro dentro il quadrato del tessuto. Distanza di un blocco = distanza del suo centro
# dalla giuntura interna piu' vicina (linee x = xmin + k T, y = ymin + k T, k = 1 ... n - 1). Fasce: [0,8) [8,16) [16,40) [40,100) [100,250].
# J-1: UMI medio per blocco per fascia. J-2: quota di trascritti catturati nati da cellule di un altro pezzo (bordo), per fascia.
# J-3: coppie di blocchi da lati opposti di una giuntura (centri a distanza < 8 um, T = 1000; < 16 um, escluso il blocco tagliato, T = 500):
# correlazione di Pearson degli UMI e coseno medio dei vettori dei geni.
# Uso: Rscript --vanilla tools/S2.6_measures.R <rundir> <T_um> [<T_um> ...]   -> <rundir>/measures_T<T>.csv, pairs_T<T>.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(Matrix))
a <- commandArgs(TRUE); dir <- a[1]; Ts <- as.numeric(a[-1])
j <- jsonlite::fromJSON(file.path(dir, "run.json")); g <- j$grid
u2 <- readBin(file.path(dir, "u2.bin"), "integer", g$n_row * g$n_col, size = 4L)
own8 <- readBin(file.path(dir, "own8.bin"), "integer", 1e8, size = 4L); dup8 <- readBin(file.path(dir, "dup8.bin"), "integer", 1e8, size = 4L)
c8 <- readRDS(file.path(dir, "c8.rds"))
nr8 <- g$n_row %/% 4L; nc8 <- g$n_col %/% 4L; n8 <- nr8 * nc8
id <- seq_len(g$n_row * g$n_col); b8 <- as.integer((((id - 1L) %/% g$n_col) %/% 4L) * nc8 + ((id - 1L) %% g$n_col) %/% 4L + 1L)
U8 <- tabulate(rep(b8, u2), n8)
i8 <- (seq_len(n8) - 1L) %/% nc8; j8 <- (seq_len(n8) - 1L) %% nc8
xc <- g$x0 + (j8 + 0.5) * 8; yc <- g$y0 + (i8 + 0.5) * 8
xmin <- g$x0 + 280; ymin <- g$y0 + 280
side <- g$n_col * 2 - 560                                   # quadrati di lato multiplo di 8 um (1000, 2000, 3000)
inside <- xc > xmin & xc < xmin + side & yc > ymin & yc < ymin + side
bands <- c(0, 8, 16, 40, 100, 250.0001)
for (Tu in Ts) {
  ks <- seq_len(round(side / Tu) - 1L); sx <- xmin + ks * Tu; sy <- ymin + ks * Tu
  dx <- vapply(xc, function(v) min(abs(v - sx)), 0); dy <- vapply(yc, function(v) min(abs(v - sy)), 0)
  d <- pmin(dx, dy); bd <- cut(d, bands, right = FALSE, labels = c("0-8", "8-16", "16-40", "40-100", "100-250"))
  k <- inside & !is.na(bd)
  m <- do.call(rbind, lapply(levels(bd), function(l) { w <- k & bd == l
    data.frame(T_um = Tu, band = l, n_blocks = sum(w), mean_umi = mean(U8[w]), sd_umi = sd(U8[w]),
               dup_frac = sum(dup8[w]) / max(1, sum(own8[w]) + sum(dup8[w]))) }))
  m <- rbind(m, data.frame(T_um = Tu, band = "all", n_blocks = sum(inside), mean_umi = mean(U8[inside]), sd_umi = sd(U8[inside]),
                           dup_frac = sum(dup8[inside]) / max(1, sum(own8[inside]) + sum(dup8[inside]))))
  m$tag <- j$tag; m$mode <- j$mode
  write.csv(m, file.path(dir, sprintf("measures_T%g.csv", Tu)), row.names = FALSE)
  # coppie attraverso le giunture verticali (x) e orizzontali (y)
  lim <- if (abs(Tu / 8 - round(Tu / 8)) < 1e-9) 8 else 16
  P <- list()
  for (s in sx) { L <- which(inside & xc < s & s - xc < lim & dy > 8); R <- L + (if (lim == 8) 1L else 2L)
    ok <- R <= n8 & inside[R]; P[[length(P) + 1]] <- cbind(L[ok], R[ok]) }
  for (s in sy) { L <- which(inside & yc < s & s - yc < lim & dx > 8); R <- L + (if (lim == 8) nc8 else 2L * nc8)
    ok <- R <= n8 & inside[R]; P[[length(P) + 1]] <- cbind(L[ok], R[ok]) }
  P <- do.call(rbind, P)
  A <- c8[, P[, 1], drop = FALSE]; B <- c8[, P[, 2], drop = FALSE]
  cosv <- Matrix::colSums(A * B) / sqrt(Matrix::colSums(A * A) * Matrix::colSums(B * B))
  pr <- data.frame(T_um = Tu, n_pairs = nrow(P), pearson_umi = cor(U8[P[, 1]], U8[P[, 2]]), mean_cosine = mean(cosv, na.rm = TRUE),
                   tag = j$tag, mode = j$mode)
  write.csv(pr, file.path(dir, sprintf("pairs_T%g.csv", Tu)), row.names = FALSE)
}
cat("misure scritte:", dir, "\n")
