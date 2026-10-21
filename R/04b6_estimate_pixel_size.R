# R/04b6_estimate_pixel_size.R - Step 1, S1.5: stima della scala in µm per pixel (design §3, C-A; BL-016, BL-029).
# Metodi: "metadata" (scalefactors_json.json di Space Ranger, controllato contro la griglia dei bin; mai il tag TIFF)
# e "nuclear_spacing" (nuclei contati sull'immagine con il laplaciano della gaussiana; scala = sqrt(N / (A_px * lambda))).
# Definizioni operative in docs/steps/S1.5.md. Dipendenze: imager, jsonlite.

.ep_mut <- function() Sys.getenv("S15_MUTANT", "none")

#' Stima la scala dell'immagine (µm per pixel)
#'
#' `"nuclear_spacing"`: scurezza d = 1 - `value` sul riquadro dei pixel del cluster dominante (fondo = 0); risposta del
#' laplaciano della gaussiana normalizzato R_s = -s^2 * lap(G_s * d), s = 2^(k/4) px da 0.5 a 32; scala caratteristica
#' s* = massimo della media di max(R_s, 0) nel cluster; nuclei = massimi locali di R_s* (finestra 2*ceil(s*) + 1) sopra la
#' soglia di Otsu della risposta positiva, dentro il cluster. Con N nuclei su A_px pixel e densita' lambda del tipo
#' dominante (nuclei/µm²): scala = sqrt(N / (A_px * lambda)). IC 95 % con bootstrap su blocchi di 64 x 64 px (B = 200,
#' seme 1). Rifiuto (nessun numero) se s* < 1 px, se s* e' all'estremo superiore della griglia o se N < 30.
#'
#' @param clust pixel con `x`, `y`, `value` (grigio in [0, 1], nuclei scuri) e la colonna dei cluster.
#' @param cell_types catalogo con `cell_type` e `density` (nuclei/mm²).
#' @param cluster_to_dominant data.frame `cluster_id`, `cell_type`; NULL con un catalogo di un tipo = tutti i pixel.
#' @param method `"nuclear_spacing"` o `"metadata"`.
#' @param scalefactors_json file `scalefactors_json.json` (solo `"metadata"`).
#' @param size_check se TRUE (predefinito) rifiuta quando il raggio dei nuclei visto (sqrt(2) * s* * scala) e' fuori da un fattore 2
#'   dal raggio atteso dal catalogo (`nucleus_to_eq_ratio`^2 * 1e6 / density); senza `nucleus_to_eq_ratio` il controllo non si fa.
#' @return lista `pixel_size_um`, `method`, `n_features`, `ci`, `sigma_px`, `refused`, `reason`, `detail`.
#' @export
estimate_pixel_size <- function(clust, cell_types, cluster_to_dominant = NULL, method = c("nuclear_spacing", "metadata"),
                                scalefactors_json = NULL, value_col = "value", cluster_col = "intensity_cluster",
                                n_boot = 200L, block_px = 64L, size_check = TRUE) {
  method <- match.arg(method)
  if (method == "metadata") return(.ep_metadata(scalefactors_json))
  if (!value_col %in% names(clust)) stop("estimate_pixel_size(): manca la colonna ", value_col, call. = FALSE)
  if (is.null(cluster_to_dominant)) {
    if (nrow(cell_types) != 1L) stop("estimate_pixel_size(): con piu' tipi serve `cluster_to_dominant`.", call. = FALSE)
    dens <- cell_types$density; keep <- rep(TRUE, nrow(clust))
  } else {
    sz <- table(clust[[cluster_col]]); big <- names(sz)[which.max(sz)]
    ty <- cluster_to_dominant$cell_type[as.character(cluster_to_dominant$cluster_id) == big]
    dens <- cell_types$density[cell_types$cell_type == ty]; keep <- as.character(clust[[cluster_col]]) == big
  }
  lam <- dens / 1e6
  x <- as.integer(round(clust$x)); y <- as.integer(round(clust$y))
  x0 <- min(x); y0 <- min(y); nx <- max(x) - x0 + 1L; ny <- max(y) - y0 + 1L
  dmat <- matrix(0, nx, ny); msk <- matrix(FALSE, nx, ny)
  dmat[cbind(x[keep] - x0 + 1L, y[keep] - y0 + 1L)] <- 1 - clust[[value_col]][keep]   # fuori dal cluster dominante = 0 (FND-13)
  msk[cbind(x[keep] - x0 + 1L, y[keep] - y0 + 1L)] <- TRUE
  A_px <- sum(msk)
  im <- imager::as.cimg(dmat)
  sig <- 2^(seq(-4, 20) / 4)
  resp <- function(s) { b <- imager::isoblur(im, s, gaussian = TRUE); h <- imager::imhessian(b, c("xx", "yy"))
    -s^2 * (as.matrix(h$xx) + as.matrix(h$yy)) }
  energy <- vapply(sig, function(s) { r <- resp(s); mean(pmax(r[msk], 0)) }, 0)
  kstar <- which.max(energy); s_star <- sig[kstar]
  refuse <- function(why) list(pixel_size_um = NA_real_, method = method, n_features = NA_integer_, ci = c(NA_real_, NA_real_),
                               sigma_px = s_star, refused = TRUE, reason = why, detail = list(sigma = sig, energy = energy))
  if (.ep_mut() != "M9") {
    if (s_star < 1) return(refuse(sprintf("nuclei non risolti (scala caratteristica %.2f px < 1 px)", s_star)))
    if (kstar == length(sig)) return(refuse("scala caratteristica all'estremo superiore della griglia"))
  }
  r <- resp(s_star); w <- 2L * ceiling(s_star) + 1L
  rmax <- as.matrix(imager::dilate_square(imager::as.cimg(r), w))
  pos <- r[msk & r > 0]
  thr <- .ep_otsu(pos)
  pk <- which(r == rmax & r > thr & msk, arr.ind = TRUE)
  N <- nrow(pk)
  if (N < 30 && .ep_mut() != "M9") return(refuse(sprintf("troppo pochi nuclei rilevati (%d < 30)", N)))
  s_hat <- sqrt(N / (A_px * lam))
  # controllo di coerenza (aggiunto in S1.5 dopo CP-1, nota D-11): il raggio dei nuclei visto dal laplaciano (sqrt(2) * s*,
  # in µm con la scala stimata) deve stare entro un fattore 2 dal raggio atteso dal catalogo (rho^2 * A_ref); altrimenti rifiuto
  ty_row <- if (is.null(cluster_to_dominant)) 1L else which(cell_types$cell_type == ty)
  rho <- if ("nucleus_to_eq_ratio" %in% names(cell_types)) cell_types$nucleus_to_eq_ratio[ty_row] else NA_real_
  r_obs <- sqrt(2) * s_star * s_hat; r_exp <- if (is.finite(rho)) sqrt(rho^2 * (1e6 / dens) / pi) else NA_real_
  if (size_check && is.finite(r_exp) && .ep_mut() != "M9" && (r_obs / r_exp < 0.5 || r_obs / r_exp > 2))
    return(c(refuse(sprintf("taglia dei nuclei incoerente con il conteggio (raggio visto %.2f um, atteso %.2f um)", r_obs, r_exp))[-8],
             list(detail = list(sigma = sig, energy = energy, A_px = A_px, threshold = thr, n = N, s_count = s_hat, r_obs = r_obs, r_exp = r_exp))))
  # bootstrap su blocchi
  bi <- (pk[, 1] - 1L) %/% block_px; bj <- (pk[, 2] - 1L) %/% block_px
  mi <- which(msk, arr.ind = TRUE); key_m <- ((mi[, 1] - 1L) %/% block_px) * 1e6 + (mi[, 2] - 1L) %/% block_px
  key_p <- bi * 1e6 + bj
  blocks <- sort(unique(key_m)); nb <- tabulate(match(key_p, blocks), length(blocks)); ab <- tabulate(match(key_m, blocks), length(blocks))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(1L)
  bs <- vapply(seq_len(n_boot), function(b) { k <- sample.int(length(blocks), replace = TRUE); sqrt(sum(nb[k]) / (sum(ab[k]) * lam)) }, 0)
  list(pixel_size_um = s_hat, method = method, n_features = N, ci = unname(stats::quantile(bs, c(0.025, 0.975))),
       sigma_px = s_star, refused = FALSE, reason = "", detail = list(sigma = sig, energy = energy, A_px = A_px, threshold = thr, r_obs = r_obs, r_exp = r_exp))
}

.ep_otsu <- function(v, nbins = 256L) {
  if (!length(v)) return(Inf)
  h <- graphics::hist(v, breaks = seq(min(v), max(v) + 1e-12, length.out = nbins + 1L), plot = FALSE)
  p <- h$counts / sum(h$counts); m <- h$mids
  w0 <- cumsum(p); mu0 <- cumsum(p * m); mt <- sum(p * m)
  sb <- (mt * w0 - mu0)^2 / (w0 * (1 - w0)); sb[!is.finite(sb)] <- 0
  m[which.max(sb)]
}

.ep_metadata <- function(f) {
  if (is.null(f) || !file.exists(f)) stop("estimate_pixel_size(): serve un file scalefactors_json esistente.", call. = FALSE)
  if (.ep_mut() == "M8") {   # mutante: tag di risoluzione del TIFF (96 dpi) invece dei metadati (BL-016)
    return(list(pixel_size_um = 25400 / 96, method = "metadata", n_features = NA_integer_, ci = c(NA_real_, NA_real_),
                sigma_px = NA_real_, refused = FALSE, reason = "", detail = list()))
  }
  j <- jsonlite::fromJSON(f)
  mpp <- j$microns_per_pixel
  if (is.null(mpp) || !is.numeric(mpp) || !(mpp > 0)) stop("estimate_pixel_size(): microns_per_pixel assente in ", f, call. = FALSE)
  # controllo contro la griglia dei bin: bin_size_um / spot_diameter_fullres (pixel per bin) = µm/px
  chk <- NA_real_
  if (!is.null(j$bin_size_um) && !is.null(j$spot_diameter_fullres)) chk <- j$bin_size_um / j$spot_diameter_fullres
  if (is.finite(chk) && abs(chk / mpp - 1) > 1e-3)
    stop(sprintf("estimate_pixel_size(): microns_per_pixel %.6g e griglia dei bin %.6g non concordano.", mpp, chk), call. = FALSE)
  list(pixel_size_um = mpp, method = "metadata", n_features = NA_integer_, ci = c(mpp, mpp), sigma_px = NA_real_,
       refused = FALSE, reason = "", detail = list(bin_check = chk, file = f))
}
