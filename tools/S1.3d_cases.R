# tools/S1.3d_cases.R — S1.3d (BL-084): casi costruiti U (u1..u5) e casi di stress CP-1. Coordinate con offset 1000 µm (ulp ~1e-13,
# come i dati reali). Ogni costruttore restituisce list(key, terr, reg, pts = NULL | data.frame, note, changed) dove `changed` asserisce (BL-083)
# che la perturbazione ha cambiato l'input. Pre-registrazione results/S1.3d/S1.3d_preregistration.md; definizioni in dev_notes.
ulp1 <- function(x) 2^(floor(log2(pmax(abs(x), 2^-1022))) - 52)   # come la revisione S1.3b: ulp(0) = 2^-1074
nextk <- function(x, k) x + k * ulp1(x)
sq <- function(x0, y0, x1, y1) st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
pg <- function(m, holes = list()) st_polygon(c(list(rbind(m, m[1, ])), lapply(holes, function(h) rbind(h, h[1, ]))))
mk <- function(key, cells, reg, pts = NULL, note = "", changed = TRUE, extra = list()) list(key = key, terr = st_sfc(cells), reg = st_sfc(reg), pts = pts, note = note, changed = changed, extra = extra)

# ---- u1: analitici (C-unit di S1.3b) -------------------------------------------------------------------------------
u1_cases <- function() {
  out <- lapply(c(1, 4, 16), function(k) mk(sprintf("U1_sq_k%02d", k), list(sq(1000, 0, 1001, 1), pg(rbind(c(1001, 0), c(1002, 0), c(1002, 1), c(nextk(1001, k), 1)))),
                                            sq(1000, 0, 1002, 1), note = "due quadrati, vertice (1001,1) del secondo spostato di k ulp", changed = nextk(1001, k) != 1001))
  out[[4]] <- mk("U1_T", list(sq(1000, 0, 1002, 2), pg(rbind(c(1002, 0), c(1004, 0), c(1004, 1), c(nextk(1002, 3), 1))), pg(rbind(c(nextk(1002, 7), 1), c(1004, 1), c(1004, 2), c(1002, 2)))),
                 sq(1000, 0, 1004, 2), note = "giunzione a T con vertici diversi di 3 e 7 ulp")
  h <- rbind(c(1003, 1003), c(1003, 1007), c(1007, 1007), c(1007, 1003))   # buco CW
  out[[5]] <- mk("U1_hole", list(pg(rbind(c(1000, 1000), c(1010, 1000), c(1010, 1010), c(1000, 1010)), list(h)), sq(1003, 1003, 1007, 1007)), sq(1000, 1000, 1010, 1010), note = "cella con buco + cella nel buco")
  out
}

# ---- u2: prove inverse della revisione S1.3b (dallo stato annodato, ±1/±4 ulp su vertice di bordo o interno) ----------
u2_cases <- function() {
  P1 <- read.csv("results/S1.3b/S1.3b_posthoc_P1_pairs_v0.csv"); out <- list()
  setv <- function(g, df) s3b_rebuild(df, class(g)[2])$geoms[[1]]
  for (q in seq_len(nrow(P1))) {
    cs <- P1[q, ]; G <- regen(roi_i(cs$archetype, cs$roi_id), cs$model, cs$rep); ti <- G$terr[[cs$i]]; tj <- G$terr[[cs$j]]
    vi <- s3b_geoms_df(list(ti)); vj <- s3b_geoms_df(list(tj)); Dm <- sqrt(outer(vi$x, vj$x, "-")^2 + outer(vi$y, vj$y, "-")^2)
    w <- which(Dm <= 1e-6, arr.ind = TRUE); dif <- w[Dm[w] > 0, , drop = FALSE]; idn <- w[Dm[w] == 0, , drop = FALSE]
    vf <- vj; vf$x[dif[, 2]] <- vi$x[dif[, 1]]; vf$y[dif[, 2]] <- vi$y[dif[, 1]]
    tag <- sprintf("%s_%s_%s_%02d", cs$archetype, cs$roi_id, cs$model, cs$rep)
    T0 <- G$terr; T0[[cs$j]] <- setv(tj, vf)
    out[[length(out) + 1]] <- mk(sprintf("U2_%s_fixed", tag), T0, G$reg[[1]], note = sprintf("stato annodato della coppia %d/%d", cs$i, cs$j), changed = !identical(vf, vj))
    for (kind in c("boundary", "interior")) { b <- if (kind == "boundary") dif[1, 2] else idn[1, 2]
      for (ax in c("x", "y")) for (s in c(-1, 1)) for (nu in c(1, 4)) {
        vv <- vf; vv[[ax]][b] <- vv[[ax]][b] + s * nu * ulp1(vv[[ax]][b]); Tk <- G$terr; Tk[[cs$j]] <- setv(tj, vv)
        out[[length(out) + 1]] <- mk(sprintf("U2_%s_%s_%s_%s%d", tag, kind, ax, if (s > 0) "p" else "m", nu), Tk, G$reg[[1]],
                                     note = sprintf("coppia %d/%d, %d ulp su vertice %s (%s)", cs$i, cs$j, nu, kind, ax), changed = !identical(vv, vf)) } }
  }
  out
}

# ---- u3: avversari per O2 (regione 40 x 40 µm a 1000) --------------------------------------------------------------
u3_cases <- function() {
  R <- sq(1000, 1000, 1040, 1040); out <- list()
  out[[1]] <- mk("U3_axis", unlist(lapply(0:3, function(a) lapply(0:3, function(b) sq(1000 + 10 * a, 1000 + 10 * b, 1010 + 10 * a, 1010 + 10 * b))), recursive = FALSE), R, note = "lati sulle linee dei bin")
  dm <- rbind(c(1020, 1010), c(1030, 1020), c(1020, 1030), c(1010, 1020))
  out[[2]] <- mk("U3_diamond", list(pg(dm), pg(rbind(c(1000, 1000), c(1020, 1000), c(1020, 1010), c(1010, 1020), c(1000, 1020))),
      pg(rbind(c(1020, 1000), c(1040, 1000), c(1040, 1020), c(1030, 1020), c(1020, 1010))), pg(rbind(c(1040, 1020), c(1040, 1040), c(1020, 1040), c(1020, 1030), c(1030, 1020))),
      pg(rbind(c(1020, 1040), c(1000, 1040), c(1000, 1020), c(1010, 1020), c(1020, 1030)))), R, note = "vertici sugli angoli dei bin")
  for (s in c(-1, 1)) { xs <- nextk(1008, s)
    out[[length(out) + 1]] <- mk(sprintf("U3_ulp_%s", if (s > 0) "p" else "m"), list(sq(1000, 1000, xs, 1040), sq(xs, 1000, 1040, 1040)), R, note = "lato a 1 ulp da x = 1008", changed = xs != 1008) }
  xb <- nextk(1008, 1); xt <- nextk(1008, -1)
  out[[length(out) + 1]] <- mk("U3_ulp_slant", list(pg(rbind(c(1000, 1000), c(xb, 1000), c(xt, 1040), c(1000, 1040))), pg(rbind(c(xb, 1000), c(1040, 1000), c(1040, 1040), c(xt, 1040)))), R, note = "lato da +1 a -1 ulp attorno a x = 1008")
  tri <- rbind(c(1009.98, 1009.99), c(1010.02, 1009.99), c(1010, 1010.02))
  out[[length(out) + 1]] <- mk("U3_tiny", list(pg(rbind(c(1000, 1000), c(1040, 1000), c(1040, 1040), c(1000, 1040)), list(tri[3:1, ])), pg(tri)), R, note = "cella < 1e-3 µm² sull'angolo di un bin")
  out[[length(out) + 1]] <- mk("U3_sliver", list(sq(1000, 1000, 1020, 1040), sq(1020, 1000, 1020 + 1e-7, 1040), sq(1020 + 1e-7, 1000, 1040, 1040)), R, note = "scheggia larga 1e-7 µm su x = 1020")
  th <- seq(0, 2 * pi, length.out = 65)[-65]; ci <- cbind(1020.3 + 15 * cos(th), 1020.7 + 15 * sin(th))
  out[[length(out) + 1]] <- mk("U3_big", list(pg(ci), pg(rbind(c(1000, 1000), c(1040, 1000), c(1040, 1040), c(1000, 1040)), list(ci[64:1, ]))), R, note = "cella su > 100 bin + complemento con buco")
  out
}

# ---- u4: vertice condiviso sul bordo della regione a cavallo del punto di arrotondamento (k + 1/2) g -------------------
u4_pair <- function(h, g = G_GRID) {           # xL < h < xR, arrotondati a nodi diversi della griglia g (asserito)
  for (k in c(2, 4, 8, 16, 32, 64)) { xL <- nextk(h, -k); xR <- nextk(h, k); if (round(xL / g) != round(xR / g)) return(c(xL, xR, k)) }
  stop("u4: nessuna coppia a cavallo trovata per ", h)
}
u4_cases <- function() {
  out <- list()
  for (q in 1:10) { h <- (1000e6 + q * 1e6 * 0.9 + 0.5) / 1e6 + 0; xx <- u4_pair(h); xt <- 1000 + q * 0.83
    out[[q]] <- mk(sprintf("U4_flat_%02d", q), list(pg(rbind(c(1000, 1000), c(xx[1], 1000), c(xt, 1010), c(1000, 1010))), pg(rbind(c(xx[2], 1000), c(1010, 1000), c(1010, 1010), c(xt, 1010)))),
                   sq(1000, 1000, 1010, 1010), note = sprintf("bordo piatto, vertici a %g ulp da (k+1/2)g", xx[3])) }
  for (q in 1:10) { h <- (1000e6 + q * 1e6 * 0.9 + 0.5) / 1e6; xx <- u4_pair(h); xt <- 1000 + q * 0.83; yl <- 1000 + 0.3 * (xx[1] - 1000); yr <- 1000 + 0.3 * (xx[2] - 1000)
    out[[10 + q]] <- mk(sprintf("U4_slope_%02d", q), list(pg(rbind(c(1000, 1000), c(xx[1], yl), c(xt, 1010), c(1000, 1010))), pg(rbind(c(xx[2], yr), c(1010, 1003), c(1010, 1010), c(xt, 1010)))),
                        pg(rbind(c(1000, 1000), c(1010, 1003), c(1010, 1010), c(1000, 1010))), note = sprintf("bordo inclinato 0.3, vertici a %g ulp da (k+1/2)g", xx[3])) }
  out
}

# ---- u5: punti per F4 (tassellazione 04b3 come riferimento) --------------------------------------------------------
u5_point_sets <- function() {
  lat <- expand.grid(a = 0:9, b = 0:9); L <- data.frame(x = 1001 + 2 * lat$a, y = 1001 + 2 * lat$b)
  notch <- pg(rbind(c(1000, 1000), c(1020, 1000), c(1020, 1020), c(1014, 1020), c(1014, 1010), c(1006, 1010), c(1006, 1020), c(1000, 1020)))
  inN <- lengths(st_intersects(st_as_sfc(lapply(seq_len(nrow(L)), function(k) st_point(c(L$x[k], L$y[k])))), st_sfc(notch))) > 0
  dia <- pg(rbind(c(1010, 1000), c(1020, 1010), c(1010, 1020), c(1000, 1010)))
  inD <- lengths(st_intersects(st_as_sfc(lapply(seq_len(nrow(L)), function(k) st_point(c(L$x[k], L$y[k])))), st_sfc(dia))) > 0
  set.seed(SEL_SEED + 50); th <- sort(runif(30, 0, 2 * pi)); rr <- 8 + 2 * runif(30); star <- cbind(1010 + rr * cos(th), 1010 + rr * sin(th))
  S <- st_sfc(pg(star)); inner <- st_sample(st_buffer(S, -1), 150); P <- st_coordinates(inner)
  near <- t(vapply(1:20, function(k) { e <- k %% 30 + 1; e2 <- e %% 30 + 1; m <- (star[e, ] + star[e2, ]) / 2; d <- c(1010, 1010) - m; m + 1e-9 * d / sqrt(sum(d^2)) }, numeric(2)))
  set.seed(SEL_SEED + 51); C <- data.frame(x = L$x + runif(100, -1e-12, 1e-12) * 0 + (runif(100) - 0.5) * 2e-12, y = L$y + (runif(100) - 0.5) * 2e-12)
  rv <- pg(rbind(c(1000, 1000), c(1012, 1000.7), c(1020, 1003.3), c(1020, 1020), c(1006, 1017.1), c(1000, 1012)))
  inR <- lengths(st_intersects(st_as_sfc(lapply(seq_len(nrow(L)), function(k) st_point(c(L$x[k], L$y[k])))), st_sfc(rv))) > 0
  list(list(key = "U5_lattice_notch", reg = notch, pts = L[inN, ], note = "reticolo dispari; incavo con lati collineari ai lati di Voronoi"),
       list(key = "U5_diamond", reg = dia, pts = L[inD, ], note = "regione che passa per i vertici di Voronoi (degeneri, 4 cocircolari)"),
       list(key = "U5_near", reg = pg(star), pts = data.frame(x = c(P[, 1], near[, 1]), y = c(P[, 2], near[, 2])), note = "20 punti a 1e-9 µm dal bordo"),
       list(key = "U5_cocirc", reg = sq(1000, 1000, 1020, 1020), pts = C, note = "reticolo perturbato di 1e-12 (vertici quasi degeneri)"),
       list(key = "U5_regvertex", reg = rv, pts = L[inR, ], note = "vertici della regione sui lati di Voronoi"))
}

# ---- CP-1: stress (60 nulli di S1.3b CP-2, stessi seed) -----------------------------------------------------------
cp1_selection <- function(ns13) { set.seed(SEL_SEED + 1); ok <- ns13[!ns13$fail, ]
  do.call(rbind, lapply(split(ok, ok$archetype), function(d) d[sample.int(nrow(d), 10), c("archetype", "roi_id", "model", "rep")])) }
job_index13b <- function(i, model, rp) (i - 1L) * 60L + (match(model, MODELS) - 1L) * 20L + rp
#' T = territori perturbati con s3b_perturb (seed di S1.3b CP-2); P = generatori perturbati U(-δ, δ) e ritassellati con 04b3
#' (un generatore che uscirebbe dalla regione non è perturbato; contato).
cp1_cases_one <- function(A, roi, model, rp) {
  i <- roi_i(A, roi); G <- regen(i, model, rp); out <- list(); base <- sprintf("%s_%s_%s_%02d", A, roi, model, rp)
  for (dl in c(1e-12, 1e-10, 1e-8)) { tag <- sprintf("d%02d", round(-log10(dl)))
    tp <- s3b_perturb(G$terr, dl, SEL_SEED + 7 * job_index13b(i, model, rp) + round(-log10(dl)))
    out[[length(out) + 1]] <- mk(sprintf("CP1T_%s_%s", base, tag), tp, G$reg[[1]], note = sprintf("territori perturbati δ = %g", dl),
                                 changed = !identical(st_as_binary(tp), st_as_binary(G$terr)))
    set.seed(SEL_SEED + 9 * job_index13b(i, model, rp) + round(-log10(dl)))
    x2 <- G$p$x + runif(length(G$p$x), -dl, dl); y2 <- G$p$y + runif(length(G$p$y), -dl, dl)
    inside <- lengths(st_intersects(st_as_sfc(lapply(seq_along(x2), function(k) st_point(c(x2[k], y2[k])))), G$reg)) > 0
    x2[!inside] <- G$p$x[!inside]; y2[!inside] <- G$p$y[!inside]
    s <- s13_tessellate(x2, y2, G$w, backend = "geos")
    out[[length(out) + 1]] <- mk(sprintf("CP1P_%s_%s", base, tag), s$tv$cell_territories, G$reg[[1]], pts = data.frame(x = x2, y = y2),
                                 note = sprintf("generatori perturbati δ = %g (%d non perturbati: uscirebbero)", dl, sum(!inside)), changed = any(x2 != G$p$x | y2 != G$p$y),
                                 extra = list(n_lost = s$df$n_lost, n_gained = s$df$n_gained))
  }
  out
}
# ---- S: pattern C-perf (4 000 x 4 000 µm, densità A4, RSA regola, seed 42, 449 536 celle) ---------------------------
s_case <- function() {
  s <- 4000; reg <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = s^2), region_polygons = st_sfc(sq(0, 0, s, s)))
  o <- seed_centroids(reg, data.frame(cell_type = "a4", density = 28096), data.frame(cluster_id = 1, cell_type = "a4", fraction = 1), random_seed = 42, verbose = FALSE)
  tv <- tessellate_voronoi(o$centroids, reg, verbose = FALSE); stopifnot(nrow(o$centroids) == 449536)
  cen <- o$centroids[order(o$centroids$cell_id), ]
  mk("S_cperf", tv$cell_territories, sq(0, 0, s, s), pts = data.frame(x = cen$x, y = cen$y), note = "pattern C-perf di S1.2/S1.3",
     extra = list(n_lost = tv$territory_df$n_pieces_lost, n_gained = tv$territory_df$n_pieces_gained, t_tessellate_04b3 = tv$info$elapsed_s))
}
