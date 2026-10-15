
# review/code/adv_cases.R — casi avversari per power_diagram() (revisore "code", S1.4b). Scrive solo in /tmp/code_adv.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow) })
L <- new.env(); sys.source("tools/S1.4b_laguerre.R", envir = L)
OD <- "/tmp/code_adv"; dir.create(OD, showWarnings = FALSE)
sq <- function(a, b, c, d) rbind(c(a, b), c(c, b), c(c, d), c(a, d), c(a, b))
cases <- list(); log <- list()
add <- function(id, x, y, w, reg) cases[[id]] <<- list(x = x, y = y, w = w, reg = reg)
set.seed(101)
# C1 pesi fortemente eterogenei (R lognormale sdlog 1.2), molti nascosti
n <- 300; R <- exp(rnorm(n, log(3), 1.2)); add("C1_hetero", runif(n, 0, 200), runif(n, 0, 200), R^2, st_polygon(list(sq(0, 0, 200, 200))))
# C1b estremo: w in [0, 1e4] uniforme
n <- 200; add("C1b_extreme", runif(n, 0, 200), runif(n, 0, 200), runif(n, 0, 1e4), st_polygon(list(sq(0, 0, 200, 200))))
# C1c pesi negativi e positivi grandi
n <- 250; add("C1c_negpos", runif(n, 0, 200), runif(n, 0, 200), rnorm(n, 0, 200), st_polygon(list(sq(0, 0, 200, 200))))
# C2 quasi coincidenti: 150 casuali + coppie a distanza 1e-2 .. 0 con pesi uguali e diversi
x <- runif(150, 0, 150); y <- runif(150, 0, 150); w <- runif(150, 0, 20)
offs <- c(1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 0)
for (k in seq_along(offs)) for (eqw in c(TRUE, FALSE)) { i <- 10 * k + eqw; x <- c(x, x[i] + offs[k]); y <- c(y, y[i]); w <- c(w, if (eqw) w[i] else w[i] + 0.5) }
add("C2_near", x, y, w, st_polygon(list(sq(0, 0, 150, 150))))
# C3a tutti collineari (n = 20 > 4)
add("C3a_allcollinear", seq(1, 96, 5), seq(1, 96, 5) / 3 + 10, runif(20, 0, 5), st_polygon(list(sq(0, 0, 100, 100))))
# C3b 40 collineari + 80 casuali, pesi casuali
add("C3b_somecollinear", c(seq(2, 98, length.out = 40), runif(80, 0, 100)), c(rep(50, 40), runif(80, 0, 100)), runif(120, 0, 10), st_polygon(list(sq(0, 0, 100, 100))))
# C3c reticolo esagonale, pesi uguali (degenere: 3 punti cocircolari per faccia, 6 per vertice) + collineari a righe
g <- expand.grid(i = 0:14, j = 0:14); add("C3c_hexlattice", 5 + 6 * g$i + 3 * (g$j %% 2), 5 + 6 * sqrt(3) / 2 * g$j, rep(4, nrow(g)), st_polygon(list(sq(0, 0, 95, 82))))
# C3d reticolo quadrato con pesi alternati (scacchiera, cocircolare con pesi)
g <- expand.grid(i = 0:11, j = 0:11); add("C3d_checker", 4 + 8 * g$i, 4 + 8 * g$j, ifelse((g$i + g$j) %% 2 == 0, 10, 2), st_polygon(list(sq(0, 0, 96, 96))))
# C4 nascosti multipli: 12 grappoli di 5 punti a peso piccolo dentro anelli di 6 punti a peso grande
x <- c(); y <- c(); w <- c()
for (cl in 1:12) { cx <- runif(1, 20, 180); cy <- runif(1, 20, 180); a <- 2 * pi * (0:5) / 6
  x <- c(x, cx + 8 * cos(a), cx + rnorm(5, 0, 1.5)); y <- c(y, cy + 8 * sin(a), cy + rnorm(5, 0, 1.5)); w <- c(w, rep(60, 6), runif(5, -30, 5)) }
x <- c(x, runif(100, 0, 200)); y <- c(y, runif(100, 0, 200)); w <- c(w, runif(100, 0, 10))
add("C4_multihidden", x, y, w, st_polygon(list(sq(0, 0, 200, 200))))
# C5 regione multipla con buchi; generatori anche fuori dalla regione
p1 <- st_polygon(list(sq(0, 0, 100, 100), sq(30, 30, 50, 50)[5:1, ], sq(60, 60, 80, 70)[5:1, ]))
p2 <- st_polygon(list(sq(120, 10, 200, 90), sq(150, 40, 170, 60)[5:1, ]))
x <- runif(260, -10, 210); y <- runif(260, -10, 110)
add("C5_multihole", x, y, runif(260, 0, 30), st_multipolygon(list(p1, p2)))
# C6 coordinate ~1e4 µm, pesi eterogenei, ritaglio
n <- 300; R <- exp(rnorm(n, log(3), 0.8)); xx <- runif(n, 0, 200); yy <- runif(n, 0, 200)
add("C6_off1e4", xx + 1e4, yy + 1e4, R^2, st_polygon(list(sq(1e4 + 5, 1e4 + 5, 1e4 + 195, 1e4 + 195))))
add("C6_off0", xx, yy, R^2, st_polygon(list(sq(5, 5, 195, 195))))
add("C6b_off1e5", xx + 1e5, yy + 1e5, R^2, st_polygon(list(sq(1e5 + 5, 1e5 + 5, 1e5 + 195, 1e5 + 195))))
# C7 celle che non contengono il generatore (peso grande accanto a piccolo)
n <- 150; x <- runif(n, 0, 120); y <- runif(n, 0, 120); w <- ifelse(runif(n) < 0.15, 400, 0)
add("C7_gen_outside", x, y, w, st_polygon(list(sq(0, 0, 120, 120))))
# C8 n piccoli (ramo n <= 4)
add("C8a_n3_hidden", c(0, 10, 5), c(0, 0, 0.1), c(50, 50, -40), st_polygon(list(sq(-20, -20, 30, 20))))
add("C8b_n4_hidden", c(0, 10, 0, 10), c(0, 0, 10, 10), c(0, 0, 0, -500), st_polygon(list(sq(-5, -5, 15, 15))))
add("C8c_n5_hidden", c(0, 10, 0, 10, 5), c(0, 0, 10, 10, 5), c(30, 30, 30, 30, -20), st_polygon(list(sq(-5, -5, 15, 15))))
add("C8d_n1", 3, 3, 1, st_polygon(list(sq(0, 0, 10, 10))))
# C9 RSA molto denso con dischi disgiunti e raggi molto diversi (teorema) — controlla contenimento
X <- numeric(0); Y <- numeric(0); Rr <- numeric(0)
for (t in 1:60000) { if (length(X) >= 400) break
  px <- runif(1, 2, 148); py <- runif(1, 2, 148); pr <- exp(runif(1, log(0.2), log(12)))
  if (!length(X) || all(sqrt((X - px)^2 + (Y - py)^2) > (Rr + pr) / 0.99)) { X <- c(X, px); Y <- c(Y, py); Rr <- c(Rr, pr) } }
add("C9_rsa_theorem", X, Y, Rr^2, st_polygon(list(sq(-15, -15, 165, 165))))
summ <- list()
for (id in names(cases)) { z <- cases[[id]]
  for (snap in c(TRUE, FALSE)) {
    t0 <- proc.time()[[3]]
    pd <- tryCatch(L$power_diagram(z$x, z$y, z$w, st_sfc(z$reg), snap = snap), error = function(e) conditionMessage(e))
    tag <- if (snap) "snap" else "exact"
    if (is.character(pd)) { summ[[length(summ) + 1]] <- data.frame(case = id, mode = tag, n = length(z$x), error = substr(pd, 1, 200), n_hidden = NA, n_multi = NA, cl3_rel = NA, c3_viol = NA); next }
    A0 <- geos_area(pd$region)
    c3 <- NA
    if (id == "C9_rsa_theorem") { nuc <- L$nucleus_polys(z$x, z$y, sqrt(z$w)); ok <- !geos_is_empty(pd$cells); c3 <- sum(!ok) + sum(!geos_within(nuc[ok], pd$cells[ok])) }
    summ[[length(summ) + 1]] <- data.frame(case = id, mode = tag, n = length(z$x), error = "", n_hidden = pd$n_hidden, n_multi = pd$n_multi,
                                           cl3_rel = abs(sum(pd$area) - A0) / A0, c3_viol = c3)
    write_parquet(data.frame(i = seq_along(z$x), x = z$x, y = z$y, w = z$w, hidden = seq_along(z$x) %in% pd$hidden, wkb = I(geos_write_wkb(pd$cells))),
                  file.path(OD, sprintf("%s_%s.parquet", id, tag)))
    writeBin(geos_write_wkb(pd$region)[[1]], file.path(OD, sprintf("%s_%s_region.wkb", id, tag)))
  }
}
s <- do.call(rbind, summ); print(s, right = FALSE); write.csv(s, file.path(OD, "R_summary.csv"), row.names = FALSE)
