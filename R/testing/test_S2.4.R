# R/testing/test_S2.4.R - S2.4: batteria di test di simulate_stripes(), capture_and_bin() e della seconda coda di
# place_and_move_transcripts() (cancello G-1, parti a-g). Pre-registrazione: docs/steps/S2.4.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S24_RSRC=<cartella con R/>] Rscript --vanilla R/testing/test_S2.4.R
# I mutanti sono copie modificate del codice (tools/S2.4_mutants.sh) passate con S24_RSRC.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
SRC <- Sys.getenv("S24_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
err <- function(expr) inherits(tryCatch({ expr; NULL }, error = function(e) e), "error")
bonf <- function(m) stats::qnorm(1 - 0.0005 / max(1, m))

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
P0 <- list(r = 2, sigma_in = 0.7, sigma_e = 4, D = 6, w = 0.1, ell = 30)
PS <- list(parity_row = 1.5, parity_col = 1.7, sd_alpha = 0.2, sd_beta = 0.25, rho_alpha = -0.3, rho_beta = -0.45, b = 0.5)

# ================= (a) conservazione
La <- mk_layer(8, 6, 10, 4, hole = c(30, 20, 50, 40)); Ha <- build_hd_grid(La, margin_um = 16)
Sa <- simulate_stripes(Ha, PS, seed = 3)
Pa <- mk_prof(La, 3000)
ma <- place_and_move_transcripts(La, Ha, Pa, c(P0, list(w2 = 0.2, ell2 = 15)), seed = 7, field = FALSE)
ca <- capture_and_bin(ma, Ha, Sa, PS, seed = 9, overgeneration = Sa$M)
X <- as.matrix(Pa$counts); cl <- ca$truth$cells; gl <- ca$truth$genes
ck("a1 per cellula: generati = profili", identical(as.numeric(cl$n_generated), as.numeric(colSums(X))))
ck("a2 per cellula: generati = catturati + scartati + persi", all(cl$n_generated == cl$n_captured + cl$n_thinned + cl$n_lost) && all(cl$n_thinned >= 0))
ck("a3 per gene: generati = profili = catturati + scartati + persi", identical(as.numeric(gl$n_generated), as.numeric(rowSums(X))) &&
     all(gl$n_generated == gl$n_captured + gl$n_thinned + gl$n_lost))
ck("a4 counts_2um per gene = catturati + fondo", identical(as.numeric(Matrix::rowSums(ca$counts_2um)), as.numeric(gl$n_captured + gl$n_background)))
nb <- Ha$grid$n_row * Ha$grid$n_col
ob <- numeric(nb); tb <- tapply(ca$truth$origin$n, ca$truth$origin$bin_id, sum); ob[as.integer(names(tb))] <- tb
bgv <- numeric(nb); bgv[ca$truth$background$bin_id] <- ca$truth$background$n
ck("a5 counts_2um per quadratino = verita' di origine + fondo", identical(as.numeric(Matrix::colSums(ca$counts_2um)), ob + bgv), sprintf("fondo %d", sum(bgv)))
ck("a6 nati nel nucleo catturati in [0, catturati]", all(ca$truth$origin$n_born_nucleus >= 0 & ca$truth$origin$n_born_nucleus <= ca$truth$origin$n))
i8 <- ((seq_len(nb) - 1) %/% Ha$grid$n_col) %/% 4; j8 <- ((seq_len(nb) - 1) %% Ha$grid$n_col) %/% 4
b8 <- i8 * (Ha$grid$n_col %/% 4) + j8 + 1
M8 <- as.matrix(ca$counts_2um %*% Matrix::sparseMatrix(i = seq_len(nb), j = b8, x = 1, dims = c(nb, max(b8))))
ck("a7 ogni quadratino di 8 um = somma dei suoi 16", identical(unname(M8), unname(as.matrix(ca$counts_8um))))
ck("a8 catturati < arrivati (le strisce assottigliano)", sum(cl$n_thinned) > 0)

# ================= (b) assottigliamento e nucleo su ingresso costruito
Hb <- Ha; nbb <- nb
fake <- function(nper, nnuc, genes = 1L) {
  tr <- data.frame(bin_id = rep(seq_len(nbb), each = genes), cell_id = 1L, gene = rep(seq_len(genes), nbb), n = as.integer(nper))
  og <- data.frame(bin_id = seq_len(nbb), cell_id = 1L, n = as.integer(nper * genes), n_born_nucleus = as.integer(nnuc))
  list(transcripts = tr, origin = og, lost = data.frame(cell_id = 1L, n_lost = 0L, n_lost_born_nucleus = 0L), genes = sprintf("g%d", seq_len(genes)),
       lost_by_gene = rep(0, genes), metadata = list(params = list()))
}
fb <- fake(200L, 60L)
cb <- capture_and_bin(fb, Hb, Sa, PS, seed = 11, overgeneration = Sa$M, background = FALSE)
pk <- Sa$c / Sa$C; kept <- numeric(nbb); kept[cb$truth$origin$bin_id] <- cb$truth$origin$n
cls <- cut(pk, stats::quantile(pk, 0:10 / 10), include.lowest = TRUE)
zb <- sapply(split(seq_len(nbb), cls), function(ix) (sum(kept[ix]) - sum(200 * pk[ix])) / sqrt(sum(200 * pk[ix] * (1 - pk[ix]))))
ck("b1 catturati = n c/C per classi di c (z, Bonferroni)", max(abs(zb)) < bonf(10), sprintf("max |z| %.2f soglia %.2f", max(abs(zb)), bonf(10)))
kn <- numeric(nbb); kn[cb$truth$origin$bin_id] <- cb$truth$origin$n_born_nucleus
zn <- (sum(kn) - sum(kept * 0.3)) / sqrt(sum(kept * 0.21 * (200 - kept) / 199))
vk <- sum((kn - kept * 0.3)^2) / sum(kept * 0.21 * (200 - kept) / 199)
ck("b2 nati nel nucleo ipergeometrici (media z e rapporto delle varianze)", abs(zn) < bonf(2) && abs(vk - 1) < 6 * sqrt(2 / nbb), sprintf("z %.2f, varianza oss/attesa %.3f", zn, vk))
set.seed(5); s0 <- .Random.seed
c0 <- capture_and_bin(fb, Hb, NULL, PS, seed = 11, overgeneration = 1, background = FALSE)
ck("b3 senza strisce e senza fondo: uscita = ingresso, generatore intatto", identical(as.numeric(Matrix::colSums(c0$counts_2um)), rep(200, nbb)) &&
     identical(c0$truth$origin$n, fb$origin$n) && identical(c0$truth$origin$n_born_nucleus, fb$origin$n_born_nucleus) && identical(.Random.seed, s0))

# ================= (c) strisce
Hc <- list(grid = list(n_row = 4000L, n_col = 4000L, bin_um = 2, x0 = 0, y0 = 0, block = 4L))
tis <- rep(c(TRUE, FALSE), each = 4000L * 2000L)                     # tessuto = prima meta' delle righe
Sc <- simulate_stripes(Hc, PS, seed = 21, tissue = tis)
i <- (seq_len(4000) - 1L); lp_r <- ifelse(i %% 2L == 1L, log(1.5), 0); lp_c <- ifelse(i %% 2L == 1L, log(1.7), 0)
Lc <- matrix(log(Sc$c * Sc$norm), 4000, 4000, byrow = TRUE)
ex <- outer(lp_r + Sc$alpha, lp_c + Sc$beta, `+`)
ck("c1 struttura esatta: log(c norm) = log p_r + log p_c + alpha + beta", max(abs(Lc - ex)) < 1e-9, sprintf("%.2e", max(abs(Lc - ex))))
ck("c2 media di c sul tessuto = 1", abs(mean(Sc$c[tis]) - 1) < 1e-12, sprintf("%.2e", mean(Sc$c[tis]) - 1))
ck("c3 M = C = max(c)", Sc$M == max(Sc$c) && Sc$C == max(Sc$c), sprintf("M %.3f", Sc$M))
ac1 <- function(x, k) { x <- x - mean(x); sum(x[-seq_len(k)] * x[seq_len(length(x) - k)]) / sum(x^2) }
okc <- TRUE; vv <- character(0)
for (rho in c(0, -0.3, -0.5)) {
  set.seed(31); xs <- replicate(100, .s24_ma1(4000L, 0.2, rho), simplify = FALSE)
  sdv <- mean(sapply(xs, stats::sd)); l1 <- mean(sapply(xs, ac1, 1)); l2 <- mean(sapply(xs, ac1, 2))
  okc <- okc && abs(sdv - 0.2) < 0.01 && abs(l1 - rho) < 0.01 && abs(l2) < 0.01
  vv <- c(vv, sprintf("rho %.2f: sd %.4f passo1 %.4f passo2 %.4f", rho, sdv, l1, l2))
}
ck("c4 media mobile di ordine 1: deviazione standard, passo 1 = rho, passo 2 = 0 (entro 0.01)", okc, paste(vv, collapse = "; "))
s <- 300 / sqrt(2) / 2; p <- as.integer(ceiling(4 * s)); d <- 0:400
kk <- exp(-0.5 * ((-p:p) / s)^2); acf_k <- sapply(c(0, 150), function(l) sum(kk[seq_len(length(kk) - l)] * kk[seq_len(length(kk) - l) + l]) / sum(kk^2))
ck("c5 campo della parita': correlazione teorica a 300 um = exp(-1/2) (entro 0.01)", abs(acf_k[2] - exp(-0.5)) < 0.01, sprintf("%.4f contro %.4f", acf_k[2], exp(-0.5)))
# varianza puntuale (un punto, molte realizzazioni): dentro una finestra piccola rispetto a 300 um la varianza spaziale e' bassa per costruzione
set.seed(41); v <- replicate(400, .s24_smooth_field(20L, 20L, 300, 2)[10, 10]); zv <- (stats::var(v) - 1) / sqrt(2 / 400)
ck("c6 campo della parita': varianza puntuale = 1 (400 realizzazioni, z)", abs(zv) < bonf(1), sprintf("varianza %.3f, z %.2f", stats::var(v), zv))
Hl <- list(grid = list(n_row = 600L, n_col = 600L, bin_um = 2), bins = list(in_tissue = rep(TRUE, 360000)))
Sl <- simulate_stripes(Hl, c(PS, list(parity_field = list(sd_row = 0.09, sd_col = 0.05, corr_um = 300))), seed = 5)
Ll <- matrix(log(Sl$c * Sl$norm), 600, 600, byrow = TRUE) - outer(Sl$alpha, Sl$beta, `+`)
odd <- (0:599) %% 2L == 1L
okl <- max(abs(Ll[!odd, !odd])) < 1e-9 && max(abs(Ll[odd, !odd] - log(1.5) - Sl$z_row[odd, !odd])) < 1e-9 && max(abs(Ll[!odd, odd] - log(1.7) - Sl$z_col[!odd, odd])) < 1e-9
ck("c7 linfonodo: il campo entra solo nelle righe e colonne dispari", okl)

# ================= (d) seconda coda
Ld <- mk_layer(20, 20, 10, 4); Hd <- build_hd_grid(Ld, margin_um = 400); Pd <- mk_prof(Ld, 400)
Pd2 <- c(P0, list(w2 = 0.25, ell2 = 12))
md <- place_and_move_transcripts(Ld, Hd, Pd, Pd2, seed = 5, field = FALSE, return_points = TRUE); pt <- md$points; nt <- nrow(pt)
f3 <- mean(pt$mode == 3L); z3 <- (f3 - 0.25) / sqrt(0.25 * 0.75 / nt)
f2 <- mean(pt$mode == 2L); z2 <- (f2 - 0.75 * 0.1) / sqrt(0.075 * 0.925 / nt)
ck("d1 frazioni dei modi: seconda coda w2, prima coda (1 - w2) w (z)", max(abs(c(z3, z2))) < bonf(2), sprintf("w2 %.4f, prima %.4f", f3, f2))
s3 <- pt$mode == 3L; rr <- sqrt((pt$x1 - pt$x0)^2 + (pt$y1 - pt$y0)^2)
ks <- suppressWarnings(stats::ks.test(rr[s3], "pgamma", shape = 2, scale = 12)$p.value)
ka <- suppressWarnings(stats::ks.test((atan2(pt$y1 - pt$y0, pt$x1 - pt$x0)[s3] + pi) / (2 * pi), "punif")$p.value)
ck("d2 raggio Gamma(2, ell2) e angolo uniforme (KS)", ks > 0.001 && ka > 0.001, sprintf("p %.3g, %.3g", ks, ka))
old <- new.env(); src_old <- tempfile(fileext = ".R")
system2("git", c("show", "c7a53de:R/05b3_place_and_move_transcripts.R"), stdout = src_old); sys.source(src_old, envir = old)
set.seed(77); r_new2 <- place_and_move_transcripts(La, Ha, Pa, c(P0, list(w2 = 0)), seed = 7, field = FALSE, frame_open = TRUE, return_points = TRUE); s_new2 <- .Random.seed
set.seed(77); r_old2 <- old$place_and_move_transcripts(La, Ha, Pa, P0, seed = 7, field = FALSE, frame_open = TRUE, return_points = TRUE); s_old2 <- .Random.seed
strip <- function(o) { o$metadata$elapsed_s <- NULL; o$metadata$params$w2 <- NULL; o }
ck("d3 w2 = 0 identico bit per bit al codice di S2.3 (c7a53de), stato del generatore compreso", identical(strip(r_new2), strip(r_old2)) && identical(s_new2, s_old2))
Lf <- mk_layer(10, 10, 10, 4); Hf <- build_hd_grid(Lf, margin_um = 16); Pf <- mk_prof(Lf, 2000)
mf <- place_and_move_transcripts(Lf, Hf, Pf, c(P0, list(w2 = 0.3, ell2 = 200)), seed = 3, field = FALSE, frame_open = TRUE, return_points = TRUE); q <- mf$points
q3 <- q$mode == 3L
inb <- all(q$x1[q3] >= 0 & q$x1[q3] <= 100 & q$y1[q3] >= 0 & q$y1[q3] <= 100)
o <- table(factor(pmin(3, floor(q$x1[q3] / 25)), 0:3), factor(pmin(3, floor(q$y1[q3] / 25)), 0:3)); chi <- stats::chisq.test(o)$p.value
ck("d4 frame_open: seconda coda mai persa, arrivi nel riquadro, uniforme dopo la riflessione (chi quadro)", sum(is.na(q$bin_id[q3])) == 0 && inb && chi > 0.001,
   sprintf("persi modo 3: %d; p %.3g", sum(is.na(q$bin_id[q3])), chi))
mg <- place_and_move_transcripts(Lf, Hf, Pf, c(P0, list(w2 = 0.3, ell2 = 200)), seed = 3, field = FALSE, frame_open = FALSE, return_points = TRUE)
lost3 <- sum(is.na(mg$points$bin_id) & mg$points$mode == 3L)
ck("d5 senza frame_open la seconda coda esce dalla griglia ed e' contata fra i persi", lost3 > 0 && sum(mg$lost$n_lost) == sum(is.na(mg$points$bin_id)), sprintf("persi modo 3: %d", lost3))

# ================= (e) attesa esatta (campo fissato, strisce, fondo, seconda coda)
Le <- mk_layer(12, 8, 10, 4, hole = c(50, 30, 70, 50)); He <- build_hd_grid(Le, margin_um = 40)
Se <- simulate_stripes(He, PS, seed = 13)
Pe <- mk_prof(Le, 600)
Pep <- list(r = 2, sigma_in = 0.7, sigma_e = 4, D = 6, w = 0, ell = 30, w2 = 0.3, ell2 = 20)
s0 <- place_and_move_transcripts(Le, He, Pe, Pep, seed = 1, field = FALSE, return_points = TRUE)
fg <- s0$field_grid; Dm <- s0$distance_matrix
set.seed(21); Fx <- s23_gamma_field(fg$n_row, fg$n_col, 0.5, 1)
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
pi_s <- ifelse(is.finite(Dm), exp(-Dm / 6), 0)
Gi <- function(u) u * stats::pnorm(u) + stats::dnorm(u)
Pax <- function(cs, b0, s, h = 0.5) { a <- cs - h / 2; outer(b0, a, function(B, A) (s * (Gi((B + 2 - A) / s) - Gi((B + 2 - A - h) / s)) - s * (Gi((B - A) / s) - Gi((B - A - h) / s))) / h) }
bx0 <- He$grid$x0 + (seq_len(He$grid$n_col) - 1) * 2; by0 <- He$grid$y0 + (seq_len(He$grid$n_row) - 1) * 2
Ecore <- matrix(0, He$grid$n_row, He$grid$n_col)
for (md_ in 1:2) { s <- c(4, 0.7)[md_]; W <- mass * if (md_ == 1) pi_s else (1 - pi_s); Ecore <- Ecore + Pax(ys, by0, s) %*% W %*% t(Pax(xs, bx0, s)) }
# seconda coda: densita' E(r) = exp(-r / ell2) / (2 pi ell2^2) al centro del quadratino ricevente per la massa di ogni sottopixel, x 4 um^2
bxc <- bx0 + 1; byc <- by0 + 1; nz <- which(mass > 0)
Et <- matrix(0, He$grid$n_row, He$grid$n_col)
for (kk in split(nz, ceiling(seq_along(nz) / 2000))) {
  sx <- gx[kk]; sy <- gy[kk]; ms <- mass[kk]
  for (r_ in seq_len(He$grid$n_row)) { dd <- sqrt(outer(sx, bxc, "-")^2 + (sy - byc[r_])^2); Et[r_, ] <- Et[r_, ] + colSums(ms * exp(-dd / 20)) / (2 * pi * 400) * 4 }
}
Cm <- matrix(Se$c, He$grid$n_row, He$grid$n_col, byrow = TRUE)
Eexp <- (0.7 * Ecore + 0.3 * Et) * Cm / Se$C + 0.5 * Cm
R_ <- 200L; obs <- matrix(0, He$grid$n_row, He$grid$n_col)
# semi delle repliche estratti a caso: place_and_move_transcripts() (S2.3) inizializza L'Ecuyer-CMRG con set.seed diretto, e semi in
# progressione aritmetica darebbero repliche legate (verifica dei flussi di S2.4)
set.seed(2024); sd_rep <- sample.int(1e8, 2 * R_)
for (k in seq_len(R_)) {
  sk <- place_and_move_transcripts(Le, He, Pe, Pep, seed = sd_rep[k], field = FALSE, field_matrix = Fx)
  ck_ <- capture_and_bin(sk, He, Se, PS, seed = sd_rep[R_ + k], overgeneration = Se$M)
  obs <- obs + matrix(Matrix::colSums(ck_$counts_2um), He$grid$n_row, He$grid$n_col, byrow = TRUE)
}
ke <- R_ * Eexp >= 20; ze <- (obs[ke] - R_ * Eexp[ke]) / sqrt(R_ * Eexp[ke])
ck("e1 media di 200 repliche = c_b [(1 - w2) regola di S2.0 + w2 E(ell2) * T] / C + b c_b (z, Bonferroni)", max(abs(ze)) < bonf(sum(ke)),
   sprintf("max |z| %.2f soglia %.2f su %d quadratini; massa della coda nella griglia %.3f", max(abs(ze)), bonf(sum(ke)), sum(ke), sum(Et) / sum(mass)))

# ================= (f) semi, processi, generatore
c1 <- capture_and_bin(ma, Ha, Sa, PS, seed = 9, overgeneration = Sa$M)
strip2 <- function(o) { o$metadata$elapsed_s <- NULL; o }
ck("f1 stesso seme stesso risultato", identical(strip2(c1), strip2(ca)))
c4a <- capture_and_bin(ma, Ha, Sa, PS, seed = 9, overgeneration = Sa$M, chunk_bins = 500L, workers = 1L)
c4b <- capture_and_bin(ma, Ha, Sa, PS, seed = 9, overgeneration = Sa$M, chunk_bins = 500L, workers = 4L)
ck("f2 1 e 4 processi identici", identical(strip2(c4a), strip2(c4b)))
set.seed(12); s0 <- .Random.seed; invisible(simulate_stripes(Ha, PS, seed = 2)); invisible(capture_and_bin(ma, Ha, Sa, PS, seed = 9, overgeneration = Sa$M))
ck("f3 generatore dell'utente ripristinato", identical(.Random.seed, s0))
fb2 <- fake(50L, 0L); half <- nbb %/% 2L
cf <- capture_and_bin(fb2, Hb, list(c = rep(1, nbb), C = 2, M = 2), PS, seed = 4, overgeneration = 2, background = FALSE, chunk_bins = half)
kf <- numeric(nbb); kf[cf$truth$origin$bin_id] <- cf$truth$origin$n
ck("f4 blocchi con lo stesso ingresso hanno estrazioni diverse", !identical(kf[1:half], kf[half + 1:half]))

empty <- list(transcripts = data.frame(bin_id = integer(0), cell_id = integer(0), gene = integer(0), n = integer(0)),
              origin = data.frame(bin_id = integer(0), cell_id = integer(0), n = integer(0), n_born_nucleus = integer(0)),
              lost = data.frame(cell_id = 1L, n_lost = 0L, n_lost_born_nucleus = 0L), genes = "g1", lost_by_gene = 1, metadata = list(params = list()))
accb <- numeric(nbb); for (k in 1:400) accb <- accb + Matrix::colSums(capture_and_bin(empty, Ha, Sa, PS, seed = 500 + k, overgeneration = Sa$M)$counts_2um)
eb <- 400 * PS$b * Sa$c; zb2 <- (accb - eb) / sqrt(eb)
ck("f5 semi consecutivi: fondo indipendente fra repliche (z per quadratino, Bonferroni)", max(abs(zb2)) < bonf(nbb), sprintf("max |z| %.2f soglia %.2f", max(abs(zb2)), bonf(nbb)))

# ================= (g) errori espliciti
ck("g1 sovra-generazione diversa da M", err(capture_and_bin(ma, Ha, Sa, PS, seed = 9, overgeneration = 1)))
ck("g2 sovra-generazione mancante", err(capture_and_bin(ma, Ha, Sa, PS, seed = 9)))
Hs <- build_hd_grid(mk_layer(3, 3, 10, 4), margin_um = 8)
ck("g3 griglia diversa", err(capture_and_bin(ma, Hs, NULL, PS, seed = 9, overgeneration = 1)))
ck("g4 |rho| > 0.5", err(simulate_stripes(Ha, modifyList(PS, list(rho_alpha = -0.6)), seed = 1)))
ck("g5 w2 fuori da [0, 1)", err(place_and_move_transcripts(La, Ha, Pa, c(P0, list(w2 = 1, ell2 = 10)), seed = 1, field = FALSE)) &&
     err(place_and_move_transcripts(La, Ha, Pa, c(P0, list(w2 = -0.1, ell2 = 10)), seed = 1, field = FALSE)))
ck("g6 ell2 mancante con w2 > 0", err(place_and_move_transcripts(La, Ha, Pa, c(P0, list(w2 = 0.1)), seed = 1, field = FALSE)))
ck("g7 seme mancante", err(simulate_stripes(Ha, PS)) && err(capture_and_bin(ma, Ha, Sa, PS, overgeneration = Sa$M)))
ck("g8 fondo acceso con b non valido", err(capture_and_bin(ma, Ha, Sa, list(b = -1), seed = 1, overgeneration = Sa$M)))

nf <- sum(!unlist(res)); cat(sprintf("\nTOTALE %d PASS, %d FAIL\n", length(res) - nf, nf))
quit(status = if (nf) 1L else 0L)
