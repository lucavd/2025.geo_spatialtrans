# review/exact/f5_cases.R — revisore esatto S1.3d, punto (a): riproduce ogni violazione O2 di F5 con polyclip installato e scrive gli interi
# che Clipper riceve (troncati come in interface.cpp: (cInt)((x - x0)/eps), e arrotondati come dice la documentazione) per i driver C++/Python.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(polyclip); library(jsonlite) })
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; O <- file.path(D, "review_exact", "f5"); dir.create(O, showWarnings = FALSE)
cat("polyclip", as.character(packageVersion("polyclip")), "\n")
# 0) il binario installato tronca o arrotonda? punto a 1.7e-6 (x/eps = 1.7): troncato -> 1e-6, arrotondato -> 2e-6
sq <- function(x0, y0, s) list(list(x = c(x0, x0 + s, x0 + s, x0), y = c(y0, y0, y0 + s, y0 + s)))
r0 <- polyclip(sq(1.7e-6, 1.7e-6, 10.2e-6), sq(0, 0, 1), op = "intersection", fillA = "nonzero", fillB = "nonzero", x0 = 0, y0 = 0, eps = 1e-6)
cat("probe: vertici di uscita min x =", min(r0[[1]]$x), " (trunc -> 1e-06, round -> 2e-06); max x =", max(r0[[1]]$x), "(trunc 1.1e-05, round 1.2e-05)\n")
cat("trunc(302/1e-6) =", sprintf("%.0f", trunc(302 / 1e-6)), " 302/1e-6 - 3.02e8 =", 302 / 1e-6 - 3.02e8, "\n")
sarea <- function(x, y) { n <- length(x); k <- c(2:n, 1); sum(x * y[k] - x[k] * y) / 2 }
pc_area <- function(R) if (!length(R)) 0 else sum(vapply(R, function(r) sarea(r$x, r$y), 0))
V <- read_parquet(file.path(D, "review_exact", "viol_all.parquet"))
V <- V[!is.na(V$k) & (grepl("F5", V$vname) | grepl("F5", V$file)) & !grepl("^o3", V$vname), ]
V$s <- ifelse(grepl("_8$", V$vname) | grepl("s8", V$file), 8, 2)
cat("violazioni F5 O2:", nrow(V), " (s2:", sum(V$s == 2), ", s8:", sum(V$s == 8), ")\n")
rings_tab <- function(key, ci) {
  t <- read_parquet(file.path(D, "geom", paste0(key, "_cells.parquet"))); t <- t[t$cell == ci, ]; t <- t[order(t$part, t$ring, t$k), ]
  out <- list()
  for (pr in unique(paste(t$part, t$ring))) { m <- t[paste(t$part, t$ring) == pr, ]; x <- m$x; y <- m$y; if (length(x) > 1 && x[1] == x[length(x)] && y[1] == y[length(y)]) { x <- x[-length(x)]; y <- y[-length(y)] }; a <- sarea(x, y)
    if ((m$ring[1] == 1 && a < 0) || (m$ring[1] > 1 && a > 0)) { x <- rev(x); y <- rev(y) }
    out[[length(out) + 1]] <- list(x = x, y = y) }
  out
}
G <- list(); INT <- list()
for (q in seq_len(nrow(V))) { r <- V[q, ]; key <- r$key; kk <- as.integer(strsplit(r$k, " ")[[1]]); ci <- kk[1]; s <- r$s
  gs <- fromJSON(file.path(D, "geom", paste0(key, "_meta.json")))$gs
  A <- rings_tab(key, ci); B <- sq(gs$X0 + s * kk[2], gs$Y0 + s * kk[3], s)
  PCx <- function(A, B, eps) polyclip(A, B, op = "intersection", fillA = "nonzero", fillB = "nonzero", x0 = 0, y0 = 0, eps = eps)
  toI <- function(R, f) lapply(R, function(z) list(x = f(z$x / 1e-6), y = f(z$y / 1e-6)))
  a_def <- pc_area(PCx(A, B, 1e-6))
  a_tr <- pc_area(PCx(toI(A, trunc), toI(B, trunc), 1)) * 1e-12
  a_rd <- pc_area(PCx(toI(A, round), toI(B, round), 1)) * 1e-12
  a_e5 <- pc_area(PCx(A, B, 1e-5)); a_e8 <- pc_area(PCx(A, B, 1e-8))
  a_auto <- pc_area(polyclip(A, B, op = "intersection", fillA = "nonzero", fillB = "nonzero"))
  G[[q]] <- data.frame(case = q, key = key, i = ci, ix = kk[2], iy = kk[3], s = s, truth = ifelse(is.na(r$area_t), 0, r$area_t), cand_reported = r$area,
                       pc_default = a_def, pc_trunc_int = a_tr, pc_round_int = a_rd, pc_eps1e5 = a_e5, pc_eps1e8 = a_e8, pc_auto = a_auto, nvert = sum(lengths(lapply(A, `[[`, "x"))))
  for (nm in c("A", "B")) { P <- if (nm == "A") A else B
    for (rr in seq_along(P)) INT[[length(INT) + 1]] <- data.frame(case = q, poly = nm, ring = rr, k = seq_along(P[[rr]]$x),
        xt = sprintf("%.0f", trunc(P[[rr]]$x / 1e-6)), yt = sprintf("%.0f", trunc(P[[rr]]$y / 1e-6)),
        xr = sprintf("%.0f", round(P[[rr]]$x / 1e-6)), yr = sprintf("%.0f", round(P[[rr]]$y / 1e-6)), x = sprintf("%.17g", P[[rr]]$x), y = sprintf("%.17g", P[[rr]]$y)) }
}
G <- do.call(rbind, G); INT <- do.call(rbind, INT)
write.csv(G, file.path(O, "f5_polyclip.csv"), row.names = FALSE); write.csv(INT, file.path(O, "f5_ints.csv"), row.names = FALSE)
cat("pc_default == cand_reported (bit):", sum(G$pc_default == G$cand_reported), " rel<=1e-9:", sum(abs(G$pc_default - G$cand_reported) <= 1e-9 * abs(G$cand_reported)), "/", nrow(G), "\n")
cat("pc_trunc_int ~ pc_default (rel 1e-9 o abs 1e-18):", sum(abs(G$pc_trunc_int - G$pc_default) <= pmax(1e-9 * abs(G$pc_default), 1e-18)), "/", nrow(G), "\n")
print(head(G[order(-abs(G$pc_default - G$truth)), ], 10))
