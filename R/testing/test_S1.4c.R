# R/testing/test_S1.4c.R — S1.4c: test di K-1/K-2 (tools/S1.4c_shape.R) e dell'operatore dei territori osservati.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S14C_MUTANT=none|MU1|MU2|MU5] Rscript --vanilla R/testing/test_S1.4c.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
source("R/04b1_extract_regions.R"); source("R/04b3_tessellate_voronoi.R"); source("tools/S1.4_mutants.R")
source("tools/S1.4c_shape.R"); L <- new.env(); sys.source("tools/S1.4b_laguerre.R", envir = L)
E <- dcg_env(); MUT <- Sys.getenv("S14C_MUTANT", "none")
res <- list(); ck <- function(name, ok, val = "") { res[[length(res) + 1]] <<- ok; cat(sprintf("%s %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, val)) }
sq <- function(x0, y0, s) geos_make_polygon(c(x0, x0 + s, x0 + s, x0), c(y0, y0, y0 + s, y0 + s))
G <- function(g) geos_set_precision(g, S14C_GRID)

# U1 nucleo che non tocca il bordo: K1 = K2 = area target (poligono a 32 vertici con fattore kp)
T1 <- G(sq(0, 0, 10)); r <- kclip_nuclei(5, 5, 50, T1, "K1", mutant = MUT)
ck("U1 K1 senza ritaglio: area = target", abs(r$df$area - 50) / 50 <= 1e-6, signif(r$df$area, 10))
r2 <- kclip_nuclei(5, 5, 50, T1, "K2", a_t = 100, phi = 1, mutant = MUT)
ck("U1 K2 senza ritaglio: area = target, 0 iterazioni", abs(r2$df$area - 50) / 50 <= 1e-6 && r2$df$n_iter == 0, r2$df$n_iter)
# U2 nucleo piu' grande del disco inscritto: K1 tronca, K2 raggiunge l'area dentro il territorio
r <- kclip_nuclei(5, 5, 90, T1, "K1", mutant = MUT)
ck("U2 K1 ritagliato: area < target", r$df$area < 90 * (1 - 1e-6), signif(r$df$area, 6))
r2 <- kclip_nuclei(5, 5, 90, T1, "K2", a_t = 100, phi = 1, mutant = MUT)
ck("U2 K2: area = target (rel <= 1e-6)", abs(r2$df$area - 90) / 90 <= 1e-6, signif(abs(r2$df$area - 90) / 90, 3))
ck("U2 K2: nucleo dentro il territorio", geos_within(r2$nuc, T1))
# U3 tetto: target 120 > A_T = 100, phi = 0.95 -> area 95, capped
r3 <- kclip_nuclei(5, 5, 120, T1, "K2", a_t = 100, phi = 0.95, mutant = MUT)
ck("U3 tetto: area = phi * A_T (territorio rientrato) dove capped", isTRUE(r3$df$capped) && abs(r3$df$area - r3$df$a_star) / r3$df$a_star <= 1e-6 && abs(r3$df$a_star - 95) < 4e-4,
   sprintf("capped=%s area=%.6f", r3$df$capped, r3$df$area))
ck("U3 C-area: area = a_star dove il tetto non e' attivo", all(abs(r3$df$area - r3$df$a_star)[!r3$df$capped] / r3$df$a_star[!r3$df$capped] <= 1e-6) &&
   all(abs(r3$df$area - r3$df$a_star) / r3$df$a_star <= 1e-6), signif(r3$df$area, 6))
# U4 due quadrati adiacenti, nuclei grandi: si toccano, non si sovrappongono; ciascuno nel proprio territorio
T2 <- G(c(sq(0, 0, 10), sq(10, 0, 10)))
r4 <- kclip_nuclei(c(5, 15), c(5, 5), c(95, 95), T2, "K2", a_t = c(100, 100), phi = 1, mutant = MUT)
ov <- geos_area(geos_intersection_prec(r4$nuc[1], r4$nuc[2], S14C_GRID))
ck("U4 C-NO: area di sovrapposizione <= 1e-9", ov <= 1e-9, signif(ov, 3))
ck("U4 C3: nuclei nei propri territori", all(geos_within(r4$nuc, T2)))
dd <- geos_distance(r4$nuc[1], r4$nuc[2])
ck("U4 contatto: distanza fra i nuclei in (0, 3e-5] (fessura 2 * inset)", dd > 0 && dd <= 3e-5, signif(dd, 3))
ck("U4 C-area", all(abs(r4$df$area - 95) / 95 <= 1e-6), signif(r4$df$area, 8))
# U5 operatore dei territori osservati (MU5): centroidi su reticolo di passo 10, generatori spostati;
# Voronoi dei centroidi = quadrati 10 x 10 -> distanza dal bordo della cella centrale = 5
reg <- list(area_um2 = 900, poly = st_sfc(st_polygon(list(rbind(c(0, 0), c(30, 0), c(30, 30), c(0, 30), c(0, 0))))))
g <- expand.grid(x = c(5, 15, 25), y = c(5, 15, 25)); set.seed(1)
gx <- g$x + runif(9, -2, 2); gy <- g$y + runif(9, -2, 2)
ob <- obs_territories(g$x, g$y, reg, E, mutant = MUT, gx = gx, gy = gy)
ck("U5 territori dei centroidi: d centrale = 5, area 100", abs(ob$d[5] - 5) < 1e-9 && abs(ob$a_t[5] - 100) < 1e-9, sprintf("d=%.4f a=%.4f", ob$d[5], ob$a_t[5]))
ck("U5 partizione: somma delle aree = 900", abs(sum(ob$a_t) - 900) < 1e-9, sum(ob$a_t))
# U6 B-8 su tre cerchi: 1-2 sovrapposti, 3 isolato
ck("U6 B-8 = 2/3", abs(b8_frac(c(0, 3, 20), c(0, 0, 0), pi * c(2, 2, 2)^2, rep(TRUE, 3), L) - 2 / 3) < 1e-12)
# U4b K1 con il disco che supera il territorio: dentro (strettamente) il proprio territorio
r4b <- kclip_nuclei(c(5, 15), c(5, 5), c(150, 150), T2, "K1", mutant = MUT)
ck("U4b C3 (K1, disco piu' grande del territorio)", all(geos_within(r4b$nuc, T2)) && all(geos_area(r4b$nuc) < 100))
r4c <- kclip_nuclei(5, 5, 100, T1, "K2", a_t = 100, phi = 1, mutant = MUT)
ck("U4c phi = 1, target = A_T: raggiunge il territorio rientrato senza errore", abs(r4c$df$area - r4c$df$a_star) / r4c$df$a_star <= 1e-6 && r4c$df$a_star < 100 && r4c$df$a_star > 100 - 4e-4,
   sprintf("a_star=%.6f area=%.6f", r4c$df$a_star, r4c$df$area))
# U7 K1/K2 su territorio non convesso (L) con generatore dentro: risultato valido, dentro il territorio
TL <- G(geos_make_polygon(c(0, 20, 20, 10, 10, 0), c(0, 0, 10, 10, 20, 20)))
r7 <- kclip_nuclei(5, 5, 250, TL, "K2", a_t = 300, phi = 1, mutant = MUT)
ck("U7 territorio non convesso: K2 valido, dentro, area = target", geos_is_valid(r7$nuc) && geos_within(r7$nuc, TL) && abs(r7$df$area - 250) / 250 <= 1e-6,
   signif(r7$df$area, 8))
nf <- sum(!unlist(res)); cat(sprintf("RESULT %d/%d PASS (mutante %s)\n", length(res) - nf, length(res), MUT))
quit(status = if (nf) 1 else 0)
