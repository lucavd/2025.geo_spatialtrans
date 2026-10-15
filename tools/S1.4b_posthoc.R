# tools/S1.4b_posthoc.R — S1.4b, analisi POST HOC (dopo i verdetti; non li cambiano).
# PH-C: le sovrapposizioni fra nuclei reali di B-8 sono un artefatto della forma circolare? Frazione di nuclei ev reali
# che intersecano un vicino, come (i) cerchio di pari area, (ii) ellisse dei momenti di pari area (ecc_n, theta_n della
# maschera; i vicini senza forma restano cerchi), (iii) ellisse ridotta del 10 % in area; per ROI, segmentatore primario.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(geos); library(parallel) })
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"; RES <- "results/S1.4b"
ARCH <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
ell <- function(x, y, A, e, th, nv = 64L) { e[!is.finite(e)] <- 0; th[!is.finite(th)] <- 0
  q <- 1 / sqrt(1 - pmin(e, 0.999)^2); a <- sqrt(A * q / pi); b <- sqrt(A / (q * pi))
  t <- 2 * pi * (seq_len(nv) - 1) / nv; lx <- outer(a, cos(t)); ly <- outer(b, sin(t))
  X <- x + lx * cos(th) - ly * sin(th); Y <- y + lx * sin(th) + ly * cos(th)
  geos_make_polygon(as.vector(t(X)), as.vector(t(Y)), feature_id = rep(seq_along(x), each = nv)) }
ovl <- function(g, ev) { tr <- geos_strtree(g); q <- geos_strtree_query(tr, g)
  i <- rep(seq_along(q), lengths(q)); j <- unlist(q); k <- i != j; i <- i[k]; j <- j[k]
  hit <- geos_intersects(g[i], g[j]); h <- rep(FALSE, length(g)); h[i[hit]] <- TRUE; mean(h[ev]) }
res <- do.call(rbind, mclapply(seq_len(nrow(rois)), function(r) { A <- rois$archetype[r]; roi <- rois$roi_id[r]; m <- ARCH[[A]]
  cc <- readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, roi, m)))$cells
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  A_all <- ifelse(cc$matched & is.finite(cc$a_n) & cc$a_n > 0, cc$a_n, cc$a_own)
  e <- ifelse(ev | (cc$matched & is.finite(cc$ecc_n)), cc$ecc_n, 0)
  data.frame(archetype = A, roi = roi, n_ev = sum(ev), ecc_med = median(cc$ecc_n[ev], na.rm = TRUE),
             circle = ovl(ell(cc$x, cc$y, A_all, rep(0, nrow(cc)), rep(0, nrow(cc))), ev),
             ellipse = ovl(ell(cc$x, cc$y, A_all, e, cc$theta_n), ev),
             ellipse_90 = ovl(ell(cc$x, cc$y, 0.9 * A_all, e, cc$theta_n), ev)) }, mc.cores = 30))
write.csv(res, file.path(RES, "S1.4b_phc.csv"), row.names = FALSE)
print(aggregate(cbind(circle, ellipse, ellipse_90, ecc_med) ~ archetype, res, mean), digits = 3)
