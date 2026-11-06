# R/testing/test_S2.3.R - S2.3: batteria di test di place_and_move_transcripts() e del campo (cancello G-1, parti a-h).
# Pre-registrazione: docs/steps/S2.3.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S23_RSRC=<cartella con R/>] Rscript --vanilla R/testing/test_S2.3.R
# I mutanti sono copie modificate del codice (tools/S2.3_mutants.sh) passate con S23_RSRC. La parte (e) chiama il prototipo di S2.0
# (tools/S2.3_test_proto.py) per lambda e la trasformata di distanza.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
SRC <- Sys.getenv("S23_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
res <- list(); NZ <- 0L
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
err <- function(expr) inherits(tryCatch({ expr; NULL }, error = function(e) e), "error")
bonf <- function(m) stats::qnorm(1 - 0.0005 / max(1, m))          # soglia di Bonferroni, errore di famiglia 0.001

# --- livelli sintetici: cellule quadrate di lato a con nucleo quadrato centrale di lato b, regione = rettangolo meno un lumo ---
mk_layer <- function(nx, ny, a, b, hole = NULL) {
  sq <- function(x0, y0, s) sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
  xy <- expand.grid(i = seq_len(nx) - 1L, j = seq_len(ny) - 1L)
  cx <- xy$i * a; cy <- xy$j * a
  if (!is.null(hole)) { k <- !(cx >= hole[1] & cx + a <= hole[3] & cy >= hole[2] & cy + a <= hole[4]); cx <- cx[k]; cy <- cy[k] }
  Tt <- sf::st_sfc(lapply(seq_along(cx), function(k) sq(cx[k], cy[k], a)))
  Nn <- sf::st_sfc(lapply(seq_along(cx), function(k) sq(cx[k] + (a - b) / 2, cy[k] + (a - b) / 2, b)))
  R <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(nx * a, 0), c(nx * a, ny * a), c(0, ny * a), c(0, 0)))))
  if (!is.null(hole)) R <- sf::st_difference(R, sf::st_sfc(sq(hole[1], hole[2], hole[3] - hole[1])))
  list(cell_df = data.frame(cell_id = 100L + seq_along(cx), cell_type = "x", stringsAsFactors = FALSE), cell_territories = Tt, cell_nuclei = Nn,
       region_polygons = R, metadata = list(pixel_size_um = 1, random_seed = 1, tissue_preset_used = "test"))
}
mk_prof <- function(layer, N, G = 5L, seed = 1) {
  set.seed(seed); n <- nrow(layer$cell_df)
  X <- vapply(seq_len(n), function(i) as.numeric(stats::rmultinom(1, N, rep(1 / G, G))), numeric(G))
  list(counts = Matrix::Matrix(X, sparse = TRUE, dimnames = list(sprintf("g%d", seq_len(G)), as.character(layer$cell_df$cell_id))))
}
P0 <- list(r = 2, sigma_in = 0.7, sigma_e = 4, D = 6, w = 0.1, ell = 30, field_s2 = 0.5, field_L = 1)

# ================= (a) conservazione
La <- mk_layer(8, 6, 10, 4, hole = c(30, 20, 50, 40)); Ha <- build_hd_grid(La, margin_um = 8); Pa <- mk_prof(La, 3000)
sa <- place_and_move_transcripts(La, Ha, Pa, P0, seed = 7, return_points = TRUE)
pt <- sa$points
X <- as.matrix(Pa$counts); Xp <- as.matrix(table(factor(pt$gene, levels = seq_len(nrow(X))), factor(pt$cell_id, levels = La$cell_df$cell_id)))
ck("a1 trascritti per cellula e gene (nati) = profili", identical(dim(Xp), dim(X)) && all(unclass(Xp) == X))
lost_p <- tapply(is.na(pt$bin_id), factor(pt$cell_id, levels = La$cell_df$cell_id), sum)
ck("a2 persi per cellula = trascritti fuori dalla griglia", identical(as.integer(lost_p), sa$lost$n_lost), sprintf("persi %d", sum(sa$lost$n_lost)))
ck("a3 nei quadratini + persi = totale della cellula", identical(as.numeric(tapply(sa$transcripts$n, factor(sa$transcripts$cell_id, levels = La$cell_df$cell_id), sum)) + sa$lost$n_lost,
                                                          as.numeric(colSums(X))))
ag <- aggregate(n ~ bin_id + cell_id + gene, data = transform(pt[!is.na(pt$bin_id), c("bin_id", "cell_id", "gene")], n = 1L), FUN = sum)
ag <- ag[order(ag$bin_id, ag$cell_id, ag$gene), ]; rownames(ag) <- NULL
ck("a4 terne = aggregazione dei trascritti", identical(ag$n, sa$transcripts$n) && identical(as.numeric(ag$bin_id), as.numeric(sa$transcripts$bin_id)) &&
     identical(as.integer(ag$gene), as.integer(sa$transcripts$gene)) && identical(as.integer(ag$cell_id), as.integer(sa$transcripts$cell_id)))
og <- aggregate(cbind(n = 1L, nb = as.integer(pt$nucleus[!is.na(pt$bin_id)])) ~ bin_id + cell_id, data = data.frame(bin_id = pt$bin_id[!is.na(pt$bin_id)], cell_id = pt$cell_id[!is.na(pt$bin_id)]), FUN = sum)
og <- og[order(og$bin_id, og$cell_id), ]
ck("a5 verita' di origine = somma delle terne, nati nel nucleo", identical(as.integer(og$n), sa$origin$n) && identical(as.integer(og$nb), sa$origin$n_born_nucleus))
ck("a6 almeno un trascritto perso (il test vede i persi)", sum(sa$lost$n_lost) > 0)

# ================= (b) nascita senza spostamento e senza campo
Lb <- mk_layer(10, 10, 10, 4); Hb <- build_hd_grid(Lb, margin_um = 8); Pb <- mk_prof(Lb, 2000)
Pz <- modifyList(P0, list(sigma_in = 1e-12, sigma_e = 1e-12, w = 0))
zmax <- 0; chi_p <- 1; inside_ok <- TRUE; flag_ok <- TRUE; mres <- 0
for (rr_ in c(0.5, 1, 2.6)) {
  sb <- place_and_move_transcripts(Lb, Hb, Pb, modifyList(Pz, list(r = rr_)), seed = 3, field = FALSE, return_points = TRUE)
  q <- sb$points; ci <- match(q$cell_id, Lb$cell_df$cell_id)
  gp <- geos::geos_make_point(q$x0, q$y0); Tg <- geos::as_geos_geometry(Lb$cell_territories)[ci]; Ng <- geos::as_geos_geometry(Lb$cell_nuclei)[ci]
  inside_ok <- inside_ok && max(geos::geos_distance(gp, Tg)) <= 1e-9
  flag_ok <- flag_ok && identical(geos::geos_intersects(gp, Ng), q$nucleus)
  mres <- max(mres, abs(q$x1 - q$x0), abs(q$y1 - q$y0))
  pn <- rr_ * 16 / (rr_ * 16 + 84); n <- nrow(q); z <- (sum(q$nucleus) - n * pn) / sqrt(n * pn * (1 - pn)); zmax <- max(zmax, abs(z))
  # uniformita': coordinate relative alla cellula in 4 x 4 quadrati di 2.5 um
  cx <- (q$x0 %% 10) %/% 2.5; cy <- (q$y0 %% 10) %/% 2.5
  ov <- function(a0, a1, b0, b1) pmax(0, pmin(a1, b1) - pmax(a0, b0))
  e <- outer(0:3, 0:3, Vectorize(function(i, j) { an <- ov(i * 2.5, i * 2.5 + 2.5, 3, 7) * ov(j * 2.5, j * 2.5 + 2.5, 3, 7); rr_ * an + (6.25 - an) }))
  o <- table(factor(cx, 0:3), factor(cy, 0:3)); ex <- e / sum(e) * n
  chi_p <- min(chi_p, stats::pchisq(sum((o - ex)^2 / ex), 15, lower.tail = FALSE))
}
ck("b1 ogni punto di nascita nel suo territorio (1e-9 um)", inside_ok)
ck("b2 nascita nel nucleo = punto nel poligono del nucleo", flag_ok)
ck("b3 senza spostamento arrivo = nascita", mres < 1e-9, sprintf("%.2e", mres))
ck("b4 quota nel nucleo = r A_n / (r A_n + A_c), r 0.5, 1, 2.6", zmax < bonf(3), sprintf("max |z| %.2f soglia %.2f", zmax, bonf(3)))
ck("b5 uniformita' in 16 celle (chi quadro, p > 0.001)", chi_p > 0.001, sprintf("p min %.3g", chi_p))

# ================= (c) spostamento su un rettangolo grande (bordo a distanze note)
Lc <- mk_layer(20, 10, 20, 2); Hc <- build_hd_grid(Lc, margin_um = 280); Pc <- mk_prof(Lc, 2000)
Pcp <- list(r = 1, sigma_in = 0.7, sigma_e = 4, D = 10, w = 0.2, ell = 30)
sc <- place_and_move_transcripts(Lc, Hc, Pc, Pcp, seed = 5, field = FALSE, return_points = TRUE)
q <- sc$points; dx <- q$x1 - q$x0; dy <- q$y1 - q$y0
n <- nrow(q); ft <- mean(q$mode == 2L)
ck("c1 frazione della coda = w", abs(ft - 0.2) / sqrt(0.2 * 0.8 / n) < bonf(8), sprintf("%.4f", ft))
rt <- sqrt(dx^2 + dy^2)[q$mode == 2L]; ks <- suppressWarnings(stats::ks.test(rt, "pgamma", shape = 2, scale = 30)$p.value)
ck("c2 raggio della coda ~ Gamma(2, ell) (KS p > 0.001)", ks > 0.001, sprintf("p %.3g", ks))
zv <- function(v, s, m) (v - s^2) / (s^2 * sqrt(2 / m))
for (md in 0:1) { s <- if (md == 1L) 4 else 0.7; k <- q$mode == md
  ck(sprintf("c3 deviazione standard per asse, modo %d (sigma %.1f)", md, s), max(abs(zv(var(dx[k]), s, sum(k))), abs(zv(var(dy[k]), s, sum(k)))) < bonf(8),
     sprintf("sd %.4f %.4f", sd(dx[k]), sd(dy[k]))) }
xc <- (floor((q$x0 - sc$field_grid$x0) / 0.5) + 0.5) * 0.5 + sc$field_grid$x0; yc <- (floor((q$y0 - sc$field_grid$y0) / 0.5) + 0.5) * 0.5 + sc$field_grid$y0
dth <- pmin(xc, 400 - xc, yc, 200 - yc) + 0.25
ck("c4 distanza dal bordo = trasformata esatta della maschera sui sottopixel", max(abs(q$d - dth)) < 1e-9, sprintf("%.2e", max(abs(q$d - dth))))
nt <- q$mode != 2L; br <- cut(q$d[nt], c(0, 2, 5, 10, 20, 40, Inf), right = FALSE)
zp <- tapply(seq_along(br), br, function(ix) { p <- mean(exp(-q$d[nt][ix] / 10)); o <- mean(q$mode[nt][ix] == 1L); (o - p) / sqrt(p * (1 - p) / length(ix) + 1e-300) })
ck("c5 frazione larga per fasce di d = exp(-d / D)", max(abs(zp[is.finite(zp)])) < bonf(8), sprintf("max |z| %.2f", max(abs(zp[is.finite(zp)]))))

# ================= (d) campo
set.seed(11); Z <- .s23_gaussian_z(1200, 1200, 1)
g <- .s23_taps(1); th <- function(k) sum(g[seq_len(length(g) - k)] * g[(k + 1):length(g)]) / sum(g^2)
cz <- function(k) mean(c(cor(as.vector(Z[, 1:(1200 - k)]), as.vector(Z[, (1 + k):1200])), cor(as.vector(Z[1:(1200 - k), ]), as.vector(Z[(1 + k):1200, ]))))
dz <- max(abs(var(as.vector(Z)) - 1), abs(cz(1) - th(1)), abs(cz(2) - th(2)), abs(cz(4) - th(4)))
ck("d1 Z: varianza 1 e correlazione a 0.5, 1, 2 um come la teorica (entro 0.02)", dz < 0.02, sprintf("%.4f", dz))
set.seed(12); Fm <- s23_gamma_field(1200, 1200, 0.5, 1)
ck("d2 F: media 1 (entro 0.02) e varianza s2 (entro il 5 %)", abs(mean(Fm) - 1) < 0.02 && abs(var(as.vector(Fm)) / 0.5 - 1) < 0.05, sprintf("media %.4f var %.4f", mean(Fm), var(as.vector(Fm))))
sub <- as.vector(Fm[seq(1, 1200, by = 20), seq(1, 1200, by = 20)])
ksF <- suppressWarnings(stats::ks.test(sub, "pgamma", shape = 2, scale = 0.5)$p.value)
ck("d3 F: margine gamma (KS su sottopixel distanti 10 um, p > 0.001)", ksF > 0.001, sprintf("p %.3g", ksF))
ck("d4 campo spento = 1", all(s23_gamma_field(10, 10, 0, 1) == 1))
tabR <- .s23_gamma_table(0.5)
pyf <- tempfile(fileext = ".py"); writeLines(c("import sys", "sys.path.insert(0, 'tools')", "import S2_3_common as S",
  "t = S.gamma_table(0.5)", "print(' '.join(repr(float(t[i])) for i in (0, 4000, 8500, 12000, 16000, 17000)))"), pyf)
pyt <- system2(".venv/bin/python", pyf, stdout = TRUE)
tp <- as.numeric(strsplit(pyt, " ")[[1]]); tr_ <- tabR[c(0, 4000, 8500, 12000, 16000, 17000) + 1]
ck("d5 tabella Z -> F identica a quella del prototipo (1e-9 relativo)", length(tp) == 6 && max(abs(tp / tr_ - 1)) < 1e-9, sprintf("%.2e", max(abs(tp / tr_ - 1))))

# ================= (e) attesa esatta e prototipo
Le <- mk_layer(12, 8, 10, 4, hole = c(50, 30, 70, 50)); He <- build_hd_grid(Le, margin_um = 40); Pe <- mk_prof(Le, 4000)
Pep <- list(r = 2, sigma_in = 0.7, sigma_e = 4, D = 6, w = 0, ell = 30)
s0 <- place_and_move_transcripts(Le, He, Pe, Pep, seed = 1, field = FALSE, return_points = TRUE)
fg <- s0$field_grid; Mt <- s0$tissue_matrix; Dm <- s0$distance_matrix
set.seed(21); Fx <- s23_gamma_field(fg$n_row, fg$n_col, 0.5, 1)
# massa per sottopixel: N_c * base * F / somma nella cellula (geometria allineata alla sottogriglia: centro del sottopixel)
xs <- fg$x0 + (seq_len(fg$n_col) - 0.5) * 0.5; ys <- fg$y0 + (seq_len(fg$n_row) - 0.5) * 0.5
gx <- matrix(xs, fg$n_row, fg$n_col, byrow = TRUE); gy <- matrix(ys, fg$n_row, fg$n_col)
ii <- floor(gx / 10); jj <- floor(gy / 10)
cid <- match(paste(ii, jj), paste(round(sapply(Le$cell_territories, function(p) min(p[[1]][, 1]))) / 10, round(sapply(Le$cell_territories, function(p) min(p[[1]][, 2]))) / 10))
cid[!(gx > 0 & gx < 120 & gy > 0 & gy < 80)] <- NA; cid <- matrix(cid, fg$n_row, fg$n_col)
innuc <- ((gx %% 10) > 3 & (gx %% 10) < 7 & (gy %% 10) > 3 & (gy %% 10) < 7)
base <- ifelse(is.na(cid), 0, ifelse(innuc, 2, 1)) * Fx
sumc <- as.numeric(tapply(as.vector(base), factor(as.vector(cid), levels = seq_len(nrow(Le$cell_df))), sum))
Nc <- colSums(as.matrix(Pe$counts))
mass <- matrix(ifelse(is.na(cid), 0, as.vector(base) * Nc[as.vector(cid)] / sumc[as.vector(cid)]), fg$n_row, fg$n_col); mass[is.na(mass)] <- 0
ck("e0 massa per sottopixel = trascritti totali", abs(sum(mass) - sum(Nc)) < 1e-6 * sum(Nc))
pi_s <- ifelse(is.finite(Dm), exp(-Dm / 6), 0)
Gi <- function(u) u * stats::pnorm(u) + stats::dnorm(u)
Pax <- function(cs, b0, s, h = 0.5) { a <- cs - h / 2; outer(b0, a, function(B, A) (s * (Gi((B + 2 - A) / s) - Gi((B + 2 - A - h) / s)) - s * (Gi((B - A) / s) - Gi((B - A - h) / s))) / h) }
bx0 <- He$grid$x0 + (seq_len(He$grid$n_col) - 1) * 2; by0 <- He$grid$y0 + (seq_len(He$grid$n_row) - 1) * 2
E <- matrix(0, He$grid$n_row, He$grid$n_col)
for (md in 1:2) { s <- c(4, 0.7)[md]; W <- mass * if (md == 1) pi_s else (1 - pi_s)
  E <- E + Pax(ys, by0, s) %*% W %*% t(Pax(xs, bx0, s)) }
R_ <- 20L; obs <- matrix(0, He$grid$n_row, He$grid$n_col)
for (k in seq_len(R_)) { sk <- place_and_move_transcripts(Le, He, Pe, Pep, seed = 100 + k, field = FALSE, field_matrix = Fx)
  tb <- tapply(sk$origin$n, sk$origin$bin_id, sum); bi <- as.integer(names(tb)); obs[cbind((bi - 1) %/% He$grid$n_col + 1, (bi - 1) %% He$grid$n_col + 1)] <- obs[cbind((bi - 1) %/% He$grid$n_col + 1, (bi - 1) %% He$grid$n_col + 1)] + tb }
k <- R_ * E >= 20; z <- (obs[k] - R_ * E[k]) / sqrt(R_ * E[k])
ck("e1 media della funzione = attesa esatta (campo fissato, z sotto Bonferroni)", max(abs(z)) < bonf(sum(k)), sprintf("max |z| %.2f soglia %.2f su %d quadratini", max(abs(z)), bonf(sum(k)), sum(k)))
# prototipo: lambda di make_blur C2 e trasformata di distanza sulla stessa maschera
td <- tempfile("s23e_"); dir.create(td)
writeBin(as.vector(t(mass)), file.path(td, "mass.bin")); writeBin(as.vector(t(1 * Mt)), file.path(td, "mask.bin"))
writeBin(as.vector(t(ifelse(is.finite(Dm), Dm, -1))), file.path(td, "dist.bin"))
bc <- expand.grid(j = seq_len(He$grid$n_col), i = seq_len(He$grid$n_row))
writeBin(c(rbind((by0[bc$i] + 1 - fg$y0), (bx0[bc$j] + 1 - fg$x0))), file.path(td, "pos.bin"))
writeLines(sprintf('{"nr": %d, "nc": %d, "sigma_e": 4, "D": 6, "w": 0, "ell": 30}', fg$n_row, fg$n_col), file.path(td, "meta.json"))
po <- system2(".venv/bin/python", c("tools/S2.3_test_proto.py", td), stdout = TRUE, stderr = TRUE)
ok_d <- any(grepl("^DIST_IDENTICAL True", po))
ck("e2 trasformata di distanza identica a quella del prototipo", ok_d, grep("^DIST", po, value = TRUE))
lam <- matrix(readBin(file.path(td, "lam.bin"), "double", n = nrow(bc)), He$grid$n_row, He$grid$n_col, byrow = TRUE)
zp <- (obs[k] - R_ * lam[k]) / sqrt(R_ * lam[k])
ck("e3 media della funzione = lambda del prototipo (z sotto Bonferroni)", max(abs(zp)) < bonf(sum(k)), sprintf("max |z| %.2f; scarto relativo massimo prototipo / esatto %.4f", max(abs(zp)), max(abs(lam[k] / E[k] - 1))))
unlink(td, recursive = TRUE)

# ================= (f) semi, processi, generatore
s1 <- place_and_move_transcripts(La, Ha, Pa, P0, seed = 9, chunk_size = 10L)
s2 <- place_and_move_transcripts(La, Ha, Pa, P0, seed = 9, chunk_size = 10L)
s4 <- place_and_move_transcripts(La, Ha, Pa, P0, seed = 9, chunk_size = 10L, workers = 4L)
s5 <- place_and_move_transcripts(La, Ha, Pa, P0, seed = 10, chunk_size = 10L)
ck("f1 stesso seme, stesso risultato", identical(s1$transcripts, s2$transcripts) && identical(s1$lost, s2$lost))
ck("f2 1 e 4 processi identici", identical(s1$transcripts, s4$transcripts) && identical(s1$origin, s4$origin))
ck("f3 seme diverso, risultato diverso", !identical(s1$transcripts, s5$transcripts))
set.seed(42); u1 <- runif(1); set.seed(42); invisible(place_and_move_transcripts(La, Ha, Pa, P0, seed = 9)); u2 <- runif(1)
ck("f4 generatore dell'utente ripristinato", identical(u1, u2) && RNGkind()[1] == "Mersenne-Twister")

# ================= (g) errori espliciti
Lg <- La; Lg$metadata$pixel_size_um <- NULL
ck("g1 scala mancante", err(place_and_move_transcripts(Lg, Ha, Pa, P0, seed = 1)))
Pg <- Pa; colnames(Pg$counts) <- rev(colnames(Pg$counts))
ck("g2 trascritti e cellule non allineati", err(place_and_move_transcripts(La, Ha, Pg, P0, seed = 1)))
bad <- list(list(sigma_e = 0), list(w = 1.5), list(D = -1), list(r = 0), list(ell = NA_real_), list(field_s2 = NA_real_), list(field_L = 0))
ck("g3 parametri non validi", all(vapply(bad, function(b) err(place_and_move_transcripts(La, Ha, Pa, modifyList(P0, b), seed = 1)), logical(1))))
Lg <- La; Lg$cell_territories[[1]] <- sf::st_polygon()
ck("g4 territorio vuoto", err(place_and_move_transcripts(Lg, Ha, Pa, P0, seed = 1)))
ck("g5 seme mancante", err(place_and_move_transcripts(La, Ha, Pa, P0)))

# ================= (h) frame_open
sh0 <- place_and_move_transcripts(Lc, Hc, Pc, Pcp, seed = 5, field = FALSE, frame_open = TRUE, return_points = TRUE)
ck("h1 frame_open: un rettangolo senza buchi non ha bordo (nessuna sfocatura larga)", all(!is.finite(sh0$points$d)) && !any(sh0$points$mode == 1L))
sh1 <- place_and_move_transcripts(La, Ha, Pa, P0, seed = 5, frame_open = TRUE, return_points = TRUE)
q <- sh1$points; xc <- (floor((q$x0 - sh1$field_grid$x0) / 0.5) + 0.5) * 0.5 + sh1$field_grid$x0; yc <- (floor((q$y0 - sh1$field_grid$y0) / 0.5) + 0.5) * 0.5 + sh1$field_grid$y0
dh <- pmax(30 - xc, xc - 50, 20 - yc, yc - 40)             # distanza (Chebyshev esterna) dal lumo [30, 50] x [20, 40]
dl <- sqrt(pmax(0, 30 - xc, xc - 50)^2 + pmax(0, 20 - yc, yc - 40)^2) + 0.25
dl2 <- ifelse(xc > 30 & xc < 50, pmax(20 - yc, yc - 40) + 0.25, ifelse(yc > 20 & yc < 40, pmax(30 - xc, xc - 50) + 0.25, NA))
ok_h <- all(abs(q$d[!is.na(dl2)] - dl2[!is.na(dl2)]) < 1e-9)
ck("h2 frame_open: la distanza viene solo dal lumo", ok_h, sprintf("d max %.1f", max(q$d)))

nP <- sum(unlist(res)); cat(sprintf("\nTOTALE: %d PASS, %d FAIL su %d\n", nP, length(res) - nP, length(res)))
