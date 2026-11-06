# R/05b4_capture_and_bin.R - Step 2, passi 5 e 6: cattura (strisce, assottigliamento, fondo) e somma a 8 um.
# Documento di progetto dello Step 2 (docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md), sezioni 7, 8, 9, 12, 13;
# sessione S2.4 (docs/steps/S2.4.md). Seconda coda dello spostamento: R/05b3 (BL-221). Dipendenze: Matrix; sum_to_8um() di R/05b1.

#' Parametri della cattura per tessuto
#'
#' Strisce: parita' dispari/pari e deviazione standard della parte fine dei fattori per riga e colonna (S2.0 PH-5, 10 iterazioni,
#' `results/S2.0/S2.0_PH5_sections.csv`, B-397...B-401); passo 1 per il generatore dalla regola R-1 di S2.4
#' (`results/S2.4/S2.4_R1_generator.csv`); campo della parita' nel linfonodo (design sezione 7). Fondo b e seconda coda (w2, ell2)
#' dalla regola R-3 di S2.4 (`results/S2.4/S2.4_R3_fits.csv`): sostituiscono il fondo del vetro lontano di R5 e i valori di DEC-051.
#' Tumore (A2) e stroma (A3) stanno sulla stessa sezione (A23).
#' @param archetype `"A1"` ... `"A6"`.
#' @return lista `parity_row`, `parity_col`, `sd_alpha`, `sd_beta`, `rho_alpha`, `rho_beta`, `parity_field` (NULL o lista `sd_row`,
#'   `sd_col`, `corr_um`), `slow` (NULL o lista `sd_alpha`, `W_alpha`, `sd_beta`, `W_beta`), `b`, `w2`, `ell2`, `section`, `archetype`.
#' @export
capture_params <- function(archetype) {
  sec <- c(A1 = "A1", A2 = "A23", A3 = "A23", A4 = "A4", A5 = "A5", A6 = "A6")
  if (!archetype %in% names(sec)) stop("capture_params(): archetipo sconosciuto: ", archetype, call. = FALSE)
  k <- sec[[archetype]]
  p <- .s24_tab[[k]]
  p$parity_field <- if (k == "A4") list(sd_row = 0.09, sd_col = 0.05, corr_um = 300) else NULL
  p$slow <- .s24_slow[[k]]
  p$section <- k; p$archetype <- archetype
  p
}

# Tabelle riempite dalle regole di S2.4 (valori con tutte le cifre dei file che li hanno prodotti; vedi capture_params()).
.s24_tab <- list(
  A1 = list(parity_row = 1.4170798953430457, parity_col = 1.750095807117609, sd_alpha = 0.1933282116659124, sd_beta = 0.1846838245230796, rho_alpha = -0.3051142894178945, rho_beta = -0.2747824688917306, b = 0.0, w2 = 0.0362255709157153, ell2 = 1825.5204420931684),
  A23 = list(parity_row = 1.577688015068367, parity_col = 1.7242280073768068, sd_alpha = 0.2209972883818394, sd_beta = 0.1908717572470308, rho_alpha = -0.3134126566799474, rho_beta = -0.2124996093518162, b = 0.1736089448280754, w2 = 0.023191356738403, ell2 = 429.4812389852601),
  A4 = list(parity_row = 0.9963845580283148, parity_col = 1.105375023614746, sd_alpha = 0.2114718928666119, sd_beta = 0.2507087192813424, rho_alpha = -0.3898006036062036, rho_beta = -0.4611807812332859, b = 0.0, w2 = 0.0363070477108733, ell2 = 1575.0424940987946),
  A5 = list(parity_row = 1.20062737537639, parity_col = 1.23141516144284, sd_alpha = 0.2419961590377538, sd_beta = 0.1803922807864428, rho_alpha = -0.3099777965083418, rho_beta = -0.296300290051404, b = 0.1474227103234672, w2 = 0.0412021259748441, ell2 = 304.59460471626346),
  A6 = list(parity_row = 0.9468670849449272, parity_col = 0.966962079996641, sd_alpha = 0.1827773294092956, sd_beta = 0.1877897009830672, rho_alpha = -0.1744003290742564, rho_beta = -0.208158682170183, b = 0.0, w2 = 0.0773522123571168, ell2 = 113.04482436271418)
)
.s24_slow <- list()

.s24_ok1 <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)

.s24_check_capture_params <- function(p) {
  for (k in c("parity_row", "parity_col", "sd_alpha", "sd_beta")) if (!.s24_ok1(p[[k]]) || p[[k]] <= 0) stop(sprintf("simulate_stripes(): `%s` deve essere un numero > 0.", k), call. = FALSE)
  for (k in c("rho_alpha", "rho_beta")) if (!.s24_ok1(p[[k]]) || abs(p[[k]]) > 0.5) stop(sprintf("simulate_stripes(): `%s` deve stare in [-0.5, 0.5] (media mobile di ordine 1).", k), call. = FALSE)
  invisible(TRUE)
}

# Media mobile di ordine 1 con deviazione standard sd e passo 1 rho: x_i = sd (e_i + theta e_{i-1}) / sqrt(1 + theta^2).
.s24_ma1 <- function(n, sd, rho) {
  th <- if (rho == 0) 0 else (1 - sqrt(max(0, 1 - 4 * rho^2))) / (2 * rho)
  e <- stats::rnorm(n + 1L)
  sd * (e[-1L] + th * e[-(n + 1L)]) / sqrt(1 + th^2)
}

# Media mobile su W valori di normali indipendenti, deviazione standard sd (parte lenta, regola R-2 di S2.4).
.s24_mavg <- function(n, sd, W) {
  W <- as.integer(W); e <- stats::rnorm(n + W - 1L)
  cs <- c(0, cumsum(e)); x <- (cs[(W + 1L):(n + W)] - cs[1:n]) / W
  sd * x * sqrt(W)
}

# Campo gaussiano sulla griglia nr x nc di passo bin_um: rumore bianco convoluto con una gaussiana di deviazione standard corr_um / sqrt(2)
# (correlazione exp(-r^2 / (2 corr_um^2))), riportato a varianza 1. Trasformata con imbottitura di 4 deviazioni standard.
.s24_smooth_field <- function(nr, nc, corr_um, bin_um) {
  s <- corr_um / sqrt(2) / bin_um; p <- as.integer(ceiling(4 * s))
  NR <- nr + 2L * p; NC <- nc + 2L * p
  dr <- pmin(0:(NR - 1L), NR - 0:(NR - 1L)); dc <- pmin(0:(NC - 1L), NC - 0:(NC - 1L))
  k <- outer(exp(-0.5 * (dr / s)^2), exp(-0.5 * (dc / s)^2))
  z <- matrix(stats::rnorm(NR * NC), NR, NC)
  f <- Re(stats::fft(stats::fft(z) * stats::fft(k), inverse = TRUE)) / (NR * NC)
  f[p + seq_len(nr), p + seq_len(nc), drop = FALSE] / sqrt(sum(k^2))
}

.s24_with_stream <- function(seed, n_streams, fun) {
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }, add = TRUE)
  # il seme si rimescola con Mersenne-Twister prima di inizializzare L'Ecuyer-CMRG: con set.seed(seme) diretto, semi in progressione
  # aritmetica danno sequenze legate linearmente (stesse estrazioni in alcune posizioni; S2.4, test e1 e verifica dei flussi)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); s1 <- sample.int(.Machine$integer.max, 1L)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(s1)
  s <- get(".Random.seed", envir = globalenv()); st <- vector("list", n_streams)
  for (k in seq_len(n_streams)) { s <- parallel::nextRNGStream(s); st[[k]] <- s }
  fun(st)
}

#' Strisce della griglia di Visium HD (fattori di cattura per quadratino)
#'
#' c_ij = p_r(i mod 2) p_c(j mod 2) exp(alpha_i + beta_j), normalizzato a media 1 sui quadratini di tessuto (design sezione 7, punto 1;
#' S2.4). alpha e beta: medie mobili di ordine 1 con la deviazione standard della parte fine e il passo 1 del generatore; parita'
#' separabile con p(pari) = 1; nel linfonodo il logaritmo del rapporto delle righe (colonne) dispari e' sommato a un campo gaussiano
#' liscio. i = riga della griglia di S2.1 (bin_id = i * n_col + j + 1).
#' @param hd uscita di `build_hd_grid()` (servono `grid` e `bins$in_tissue`).
#' @param params lista come `capture_params()`.
#' @param seed intero.
#' @param tissue vettore logico per quadratino (ordine di `bin_id`) su cui normalizzare; predefinito `hd$bins$in_tissue`.
#' @return lista `c` (fattore per quadratino, ordine di `bin_id`), `C` (= max(c)), `M` (fattore di sovra-generazione = C), `alpha`, `beta`,
#'   `slow_alpha`, `slow_beta`, `z_row`, `z_col` (campi della parita' o NULL), `norm` (media di c prima della normalizzazione), `seed`, `params`.
#' @export
simulate_stripes <- function(hd, params, seed, tissue = hd$bins$in_tissue) {
  if (missing(seed) || length(seed) != 1L || !is.finite(seed)) stop("simulate_stripes(): serve un `seed`.", call. = FALSE)
  .s24_check_capture_params(params)
  g <- hd$grid; nr <- g$n_row; nc <- g$n_col
  if (length(tissue) != nr * nc) stop("simulate_stripes(): `tissue` deve avere un valore per quadratino.", call. = FALSE)
  if (!any(tissue)) stop("simulate_stripes(): nessun quadratino di tessuto su cui normalizzare.", call. = FALSE)
  .s24_with_stream(seed, 6L, function(st) {
    assign(".Random.seed", st[[1]], envir = globalenv()); al <- .s24_ma1(nr, params$sd_alpha, params$rho_alpha)
    assign(".Random.seed", st[[2]], envir = globalenv()); be <- .s24_ma1(nc, params$sd_beta, params$rho_beta)
    sa <- sb <- NULL
    if (!is.null(params$slow)) {
      if (!is.null(params$slow$sd_alpha)) { assign(".Random.seed", st[[5]], envir = globalenv()); sa <- .s24_mavg(nr, params$slow$sd_alpha, params$slow$W_alpha); al <- al + sa }
      if (!is.null(params$slow$sd_beta)) { assign(".Random.seed", st[[6]], envir = globalenv()); sb <- .s24_mavg(nc, params$slow$sd_beta, params$slow$W_beta); be <- be + sb }
    }
    odd_r <- (seq_len(nr) - 1L) %% 2L == 1L; odd_c <- (seq_len(nc) - 1L) %% 2L == 1L
    lr <- matrix(ifelse(odd_r, log(params$parity_row), 0), nr, nc); lc <- matrix(ifelse(odd_c, log(params$parity_col), 0), nr, nc, byrow = TRUE)
    zr <- zc <- NULL
    if (!is.null(params$parity_field)) {
      pf <- params$parity_field
      assign(".Random.seed", st[[3]], envir = globalenv()); zr <- pf$sd_row * .s24_smooth_field(nr, nc, pf$corr_um, g$bin_um)
      assign(".Random.seed", st[[4]], envir = globalenv()); zc <- pf$sd_col * .s24_smooth_field(nr, nc, pf$corr_um, g$bin_um)
      lr[odd_r, ] <- lr[odd_r, ] + zr[odd_r, ]; lc[, odd_c] <- lc[, odd_c] + zc[, odd_c]
    }
    Lm <- lr + lc + outer(al, be, `+`)
    craw <- as.vector(t(exp(Lm)))                                   # ordine di bin_id: riga per riga
    nrm <- mean(craw[tissue]); cc <- craw / nrm
    list(c = cc, C = max(cc), M = max(cc), alpha = al, beta = be, slow_alpha = sa, slow_beta = sb, z_row = zr, z_col = zc,
         norm = nrm, seed = seed, params = params)
  })
}

#' Cattura dei trascritti e somma a 8 um
#'
#' Assottigliamento binomiale esatto: ogni terna (quadratino b, cellula, gene) con n trascritti arrivati diventa Binomiale(n, c_b / C);
#' i profili devono essere stati generati con il fattore di profondita' moltiplicato per M = C (`overgeneration`), cosi' l'atteso
#' catturato e' quello della scala di R4 (design sezione 7, punto 3). Per coppia (quadratino, cellula) i nati nel nucleo fra i catturati
#' sono ipergeometrici. Fondo: Poisson(b c_b) per quadratino, geni con la composizione dei trascritti generati nell'area.
#' @param moved uscita di `place_and_move_transcripts()`.
#' @param hd uscita di `build_hd_grid()` sullo stesso livello.
#' @param stripes uscita di `simulate_stripes()`, o NULL (nessuna striscia: c = 1).
#' @param params lista con `b` (fondo, UMI per quadratino; ignorato con `background = FALSE`).
#' @param seed intero.
#' @param overgeneration fattore usato per i profili; deve valere `stripes$M` (1 senza strisce).
#' @param background `TRUE` per il fondo.
#' @param chunk_bins quadratini per blocco dell'assottigliamento (fa parte della definizione del risultato).
#' @param workers processi (`parallel::mclapply`).
#' @return lista `counts_2um` (dgCMatrix geni x quadratini, ordine di bin_id), `counts_8um`, `bins` (con `c`), `truth` (`origin`,
#'   `background`, `cells`, `genes`, `layers`), `metadata`.
#' @export
capture_and_bin <- function(moved, hd, stripes, params, seed, overgeneration, background = TRUE, chunk_bins = 50000L, workers = 1L) {
  t0 <- proc.time()[["elapsed"]]
  if (missing(seed) || length(seed) != 1L || !is.finite(seed)) stop("capture_and_bin(): serve un `seed`.", call. = FALSE)
  g <- hd$grid; nb <- g$n_row * g$n_col
  tr <- moved$transcripts; og <- moved$origin; G <- length(moved$genes)
  if (is.null(tr) || is.null(og) || G == 0L) stop("capture_and_bin(): `moved` non e' un'uscita di place_and_move_transcripts().", call. = FALSE)
  if (nrow(tr) && (max(tr$bin_id) > nb || min(tr$bin_id) < 1)) stop("capture_and_bin(): quadratini fuori dalla griglia di `hd` (griglia diversa?).", call. = FALSE)
  cvec <- if (is.null(stripes)) rep(1, nb) else stripes$c
  if (length(cvec) != nb) stop("capture_and_bin(): `stripes$c` non ha un valore per quadratino di `hd`.", call. = FALSE)
  C <- if (is.null(stripes)) 1 else stripes$C
  if (missing(overgeneration) || !.s24_ok1(overgeneration) || abs(overgeneration - C) > 1e-12 * max(1, C))
    stop(sprintf("capture_and_bin(): `overgeneration` (%s) deve essere uguale a stripes$M (%.15g): i profili vanno generati con il fattore di profondita' per M.",
                 if (missing(overgeneration)) "mancante" else format(overgeneration, digits = 15), C), call. = FALSE)
  if (isTRUE(background) && (!.s24_ok1(params$b) || params$b < 0)) stop("capture_and_bin(): fondo acceso ma `b` non valido.", call. = FALSE)
  p <- cvec / C
  thin <- !is.null(stripes)
  res <- .s24_with_stream(seed, 2L, function(st) {
    # (1) assottigliamento per blocchi di quadratini (un flusso per blocco: indipendente dal numero di processi)
    kept <- tr$n; kept_nuc <- og$n_born_nucleus
    if (thin && nrow(tr)) {
      blk_t <- (tr$bin_id - 1L) %/% chunk_bins; blk_o <- (og$bin_id - 1L) %/% chunk_bins
      ub <- sort(unique(c(blk_t, blk_o)))
      # sotto-flussi del flusso 1 (prima: flussi successivi, il secondo blocco coincideva con il flusso del fondo)
      ss <- vector("list", length(ub)); s <- st[[1]]
      for (k in seq_along(ub)) { ss[[k]] <- s; s <- parallel::nextRNGSubStream(s) }
      one <- function(k) {
        assign(".Random.seed", ss[[k]], envir = globalenv())
        it <- which(blk_t == ub[k]); io <- which(blk_o == ub[k])
        kt <- stats::rbinom(length(it), tr$n[it], p[tr$bin_id[it]])
        # catturati per (quadratino, cellula) nel blocco: terne e verita' di origine hanno lo stesso ordine (quadratino, cellula),
        # quindi i gruppi consecutivi delle terne sono le righe dell'origine; poi nati nel nucleo ipergeometrici
        bt <- tr$bin_id[it]; ct <- tr$cell_id[it]; nt <- length(it)
        grp <- if (nt) cumsum(c(TRUE, bt[-1L] != bt[-nt] | ct[-1L] != ct[-nt])) else integer(0)
        first <- which(!duplicated(grp))
        if (length(first) != length(io) || any(bt[first] != og$bin_id[io]) || any(ct[first] != og$cell_id[io]))
          stop("capture_and_bin(): terne e verita' di origine non allineate.", call. = FALSE)
        kk <- as.numeric(rowsum(kt, grp, reorder = FALSE))
        kn <- stats::rhyper(length(io), og$n_born_nucleus[io], og$n[io] - og$n_born_nucleus[io], kk)
        list(it = it, kt = kt, io = io, kk = kk, kn = kn)
      }
      out <- if (workers > 1L) parallel::mclapply(seq_along(ub), one, mc.cores = workers, mc.preschedule = FALSE) else lapply(seq_along(ub), one)
      bad <- vapply(out, function(o) inherits(o, "try-error"), logical(1))
      if (any(bad)) stop("capture_and_bin(): errore in un blocco: ", as.character(out[[which(bad)[1]]]), call. = FALSE)
      kept_o <- og$n
      for (o in out) { kept[o$it] <- o$kt; kept_o[o$io] <- o$kk; kept_nuc[o$io] <- o$kn }
    } else kept_o <- og$n
    # (2) fondo
    bg <- NULL
    if (isTRUE(background)) {
      assign(".Random.seed", st[[2]], envir = globalenv())
      K <- stats::rpois(nb, params$b * cvec)
      gen_g <- as.numeric(tapply(tr$n, factor(tr$gene, levels = seq_len(G)), sum)); gen_g[is.na(gen_g)] <- 0
      gen_g <- gen_g + as.numeric(moved$lost_by_gene)
      if (sum(K) > 0) {
        gb <- sample.int(G, sum(K), replace = TRUE, prob = gen_g / sum(gen_g))
        bb <- rep(seq_len(nb), K)
        bg <- list(K = K, bin = bb, gene = gb)
      } else bg <- list(K = K, bin = integer(0), gene = integer(0))
    }
    list(kept = kept, kept_o = kept_o, kept_nuc = kept_nuc, bg = bg)
  })
  kept <- res$kept; kt <- kept > 0
  ii <- c(tr$gene[kt], if (!is.null(res$bg)) res$bg$gene); jj <- c(tr$bin_id[kt], if (!is.null(res$bg)) res$bg$bin)
  xx <- c(kept[kt], if (!is.null(res$bg)) rep(1, length(res$bg$gene)))
  X2 <- Matrix::sparseMatrix(i = as.integer(ii), j = as.integer(jj), x = as.numeric(xx), dims = c(G, nb), dimnames = list(moved$genes, NULL))
  X8 <- sum_to_8um(X2, g)
  ko <- res$kept_o > 0
  origin <- data.frame(bin_id = og$bin_id[ko], cell_id = og$cell_id[ko], n = as.integer(res$kept_o[ko]), n_born_nucleus = as.integer(res$kept_nuc[ko]))
  bgdf <- if (!is.null(res$bg)) { w <- which(res$bg$K > 0); data.frame(bin_id = w, n = as.integer(res$bg$K[w])) } else data.frame(bin_id = integer(0), n = integer(0))
  cid <- moved$lost$cell_id
  fc <- factor(tr$cell_id, levels = cid)
  arr <- as.numeric(tapply(tr$n, fc, sum)); arr[is.na(arr)] <- 0
  cap <- as.numeric(tapply(kept, fc, sum)); cap[is.na(cap)] <- 0
  cells <- data.frame(cell_id = cid, n_generated = arr + moved$lost$n_lost, n_lost = moved$lost$n_lost, n_thinned = arr - cap, n_captured = cap)
  fg <- factor(tr$gene, levels = seq_len(G))
  garr <- as.numeric(tapply(tr$n, fg, sum)); garr[is.na(garr)] <- 0; gcap <- as.numeric(tapply(kept, fg, sum)); gcap[is.na(gcap)] <- 0
  genes <- data.frame(gene = moved$genes, n_generated = garr + as.numeric(moved$lost_by_gene), n_lost = as.numeric(moved$lost_by_gene),
                      n_thinned = garr - gcap, n_captured = gcap,
                      n_background = if (!is.null(res$bg)) tabulate(res$bg$gene, G) else 0)
  bins <- hd$bins; bins$c <- cvec
  layers <- list(stripes = if (is.null(stripes)) NULL else stripes[c("alpha", "beta", "slow_alpha", "slow_beta", "z_row", "z_col", "norm", "C", "M", "seed", "params")],
                 background = list(on = isTRUE(background), b = if (isTRUE(background)) params$b else 0),
                 second_tail = moved$metadata$params[c("w2", "ell2")], overgeneration = overgeneration)
  list(counts_2um = X2, counts_8um = X8, bins = bins,
       truth = list(origin = origin, background = bgdf, cells = cells, genes = genes, layers = layers),
       metadata = list(step = "S2.4", seed = seed, background = isTRUE(background), stripes = !is.null(stripes), chunk_bins = chunk_bins,
                       n_arrived = sum(tr$n), n_captured = sum(kept), n_background = sum(bgdf$n), n_lost = sum(moved$lost$n_lost),
                       package_version = "geo.spatialtrans step2-visium-grid (S2.4)", elapsed_s = proc.time()[["elapsed"]] - t0))
}
