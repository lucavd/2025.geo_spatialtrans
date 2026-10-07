# Revisione avversariale S1.1 — 10_fixcheck.R
# Rilancia i casi avversari di 09_adversarial.R contro la versione COMMITTATA (e7f2562, copia in
# review/old_04b1_extract_regions.R) e la working copy corrente; cattura warning/errori; verifica regressioni
# (output identici) su input pieni a stride 1, mappe reali, nulli, sintetici, avversari del test.
.libPaths("/home/user/2025.geo_spatialtrans/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel) })
REV <- "/mnt/micron/geo_spatialtrans/S1.1/review"; IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
OLD <- new.env(); sys.source(file.path(REV, "old_04b1_extract_regions.R"), envir = OLD)
NEW <- new.env(); sys.source("/home/user/2025.geo_spatialtrans/R/04b1_extract_regions.R", envir = NEW)
run <- function(env, ...) {
  w <- character(0)
  v <- tryCatch(withCallingHandlers(env$extract_regions(..., verbose = FALSE),
                warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") }),
                error = function(e) structure(conditionMessage(e), class = "err"))
  list(v = v, w = w)
}
wtag <- function(w) if (!length(w)) "nessun warning" else paste0("WARNING: ", paste(sub("^extract_regions\\(\\): ", "", w), collapse = " || "))
desc <- function(z) {
  if (inherits(z$v, "err")) return(paste("ERRORE:", sub("^extract_regions\\(\\): ", "", z$v)))
  r <- z$v
  sprintf("stride %d%s, %d cluster (id %s), %d comp, %d reg, aree %s; %s", r$info$stride,
          if (!is.null(r$info$stride_source)) paste0(" [", r$info$stride_source, "]") else "",
          length(r$info$cluster_ids), paste(head(r$info$cluster_ids, 6), collapse = ","), r$info$n_components, r$info$n_regions,
          paste(head(signif(r$region_df$area_um2, 6), 4), collapse = "+"), wtag(z$w))
}
nz <- function(df) { if (is.data.frame(df) && is.numeric(df$cluster_id)) df$cluster_id <- as.numeric(df$cluster_id); df }
same_df <- function(a, b) identical(nz(a), nz(b))   # stessi valori; il tipo intero/double di cluster_id e' riportato a parte
ROWS <- list()
rec <- function(id, caso, before, after, esito, nota = "") {
  ROWS[[length(ROWS) + 1L]] <<- data.frame(id = id, caso = caso, prima = before, dopo = after, esito = esito, nota = nota)
  cat(sprintf("[%s] %s | %s\n   prima: %s\n   dopo : %s\n", esito, id, caso, before, after))
}
both <- function(...) list(o = run(OLD, ...), n = run(NEW, ...))
has_w <- function(z, pat) any(grepl(pat, z$w))
area_tot <- function(z) if (inherits(z$v, "err")) NA else sum(z$v$region_df$area_um2)

# ---- RA-17: stride ------------------------------------------------------------------------------
cl <- data.frame(x = c(1L, 4L), y = c(1L, 1L), k = 1L)
b <- both(cl, 8, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-01", "RA-17 A1: bin (1,1),(4,1) a 8 µm, stride non dichiarato", desc(b$o), desc(b$n),
    if (has_w(b$n, "stride 3 stimato") && !has_w(b$o, "stride")) "warning atteso presente" else "NON atteso",
    "Il risultato resta stride 3 / 1 componente da 1152 µm² (scelta di progetto: warning, non errore); info$stride_source = 'estimated'")
b <- both(cl, 8, cluster_col = "k")
rec("FX-02", "RA-17 A1b: stesso input, parametri d'uso (min 100, tol 0.5)", desc(b$o), desc(b$n),
    if (has_w(b$n, "stride 3 stimato") && !has_w(b$n, "simplify")) "warning atteso presente" else "NON atteso",
    "Senza stride dichiarato la regione spuria da 1152 µm² viene ancora tenuta")
b <- both(cl, 8, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k", stride = 1)
rec("FX-03", "RA-17 A1c: stesso input con stride = 1 dichiarato", desc(b$o), desc(b$n),
    if (!length(b$n$w) && b$n$v$info$n_components == 2 && identical(b$n$v$info$stride_source, "user")) "corretto, nessun warning" else "NON atteso")
cl <- data.frame(x = c(1L, 3L, 5L), y = 1L, k = 1L)
b <- both(cl, 1, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-04", "RA-17 A2: riga x = 1,3,5 a 1 µm, stride non dichiarato", desc(b$o), desc(b$n),
    if (has_w(b$n, "stride 2 stimato")) "warning atteso presente" else "NON atteso")
b2 <- both(cl, 1, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k", stride = 1)
rec("FX-05", "RA-17 A2 con stride = 1", desc(b2$o), desc(b2$n),
    if (!length(b2$n$w) && b2$n$v$info$n_components == 3) "corretto, nessun warning" else "NON atteso")
cl <- data.frame(x = c(1L, 5001L), y = c(1L, 5001L), k = 1L)
b <- both(cl, 8, simplify_tol_um = 0, cluster_col = "k")
rec("FX-06", "RA-17 A3: bin (1,1),(5001,5001) a 8 µm, stride non dichiarato", desc(b$o), desc(b$n),
    if (has_w(b$n, "stride 5000 stimato")) "warning atteso presente" else "NON atteso",
    sprintf("area totale ancora %.3g µm² senza stride dichiarato", area_tot(b$n)))
b <- both(cl, 8, simplify_tol_um = 0, cluster_col = "k", stride = 1)
rec("FX-07", "RA-17 A3 con stride = 1", desc(b$o), desc(b$n),
    if (!length(b$n$w) && b$n$v$info$n_regions == 0 && b$n$v$info$n_components == 2) "corretto, nessun warning" else "NON atteso")
cl <- data.frame(x = 3L, y = 5L, k = 1L)
b <- both(cl, 1, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-08", "RA-17 pixel singolo (adv_pixel35), stride non dichiarato", desc(b$o), desc(b$n),
    if (!length(b$n$w) && identical(b$n$v$info$stride_source, "default")) "nessun warning (corretto)" else "NON atteso")
# ---- RA-18: etichette -----------------------------------------------------------------------------
fm <- function(m) { i <- which(m != "", arr.ind = TRUE); data.frame(x = i[, 2], y = i[, 1], k = m[i], stringsAsFactors = FALSE) }
m <- matrix(rep(c("1", "1", "01", "01"), each = 2), 2, 4)
b <- both(fm(m), 10, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-09", "RA-18 A4: etichette carattere '1' e '01'", desc(b$o), desc(b$n),
    if (!inherits(b$n$v, "err") && b$n$v$info$n_regions == 2) "corretto" else "NON atteso")
cl <- fm(m); cl$k <- factor(ifelse(cl$k == "1", "1", "1.0"))
b <- both(cl, 10, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-10", "RA-18 A4b: fattore con livelli '1' e '1.0'", desc(b$o), desc(b$n),
    if (!inherits(b$n$v, "err") && b$n$v$info$n_regions == 2) "corretto" else "NON atteso")
m <- matrix(rep(c("2", "10", "1", "3"), each = 2), 2, 4)
b <- both(fm(m), 10, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
same <- same_df(b$o$v$region_df, b$n$v$region_df) && isTRUE(all.equal(as.numeric(b$o$v$info$cluster_ids), as.numeric(b$n$v$info$cluster_ids)))
rec("FX-11", "RA-18 non regressione: carattere '2','10','1','3' (conversione iniettiva)", desc(b$o), desc(b$n),
    if (same && is.numeric(b$n$v$info$cluster_ids)) "invariato (id numerici, ordine numerico)" else "NON atteso")
cl <- fm(m); cl$k <- factor(cl$k, levels = c("3", "1", "2", "10"))
b <- both(cl, 10, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
same <- same_df(b$o$v$region_df, b$n$v$region_df)
rec("FX-12", "RA-18 non regressione: fattore con livelli numerici non ordinati", desc(b$o), desc(b$n),
    if (same) "invariato" else "NON atteso")
cl <- fm(matrix(rep(c("a", "b", "a10", "b"), each = 2), 2, 4))
b <- both(cl, 10, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-13", "RA-18 non regressione: etichette testuali", desc(b$o), desc(b$n),
    if (same_df(b$o$v$region_df, b$n$v$region_df)) "invariato" else "NON atteso")
# ---- RA-19: casi limite --------------------------------------------------------------------------
cl <- data.frame(x = rep(1:3, each = 2), y = rep(1:2, 3), k = c(0L, 0L, -1L, -1L, 2L, 2L))
b <- both(cl, 5, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-14", "RA-19 cluster 0, -1, 2", desc(b$o), desc(b$n),
    if (same_df(b$o$v$region_df, b$n$v$region_df) && !length(b$n$w)) "invariato, nessun warning" else "NON atteso")
cl <- expand.grid(x = -5:4, y = -3:6); cl$k <- 1L
b <- both(cl, 2, cluster_col = "k")
rec("FX-15", "RA-19 coordinate negative (px 2, parametri d'uso)", desc(b$o), desc(b$n),
    if (same_df(b$o$v$region_df, b$n$v$region_df) && identical(st_as_binary(b$o$v$region_polygons), st_as_binary(b$n$v$region_polygons)) && !length(b$n$w)) "invariato, nessun warning" else "NON atteso")
cl <- expand.grid(x = seq(1L, 19L, 2L), y = seq(1L, 28L, 3L)); cl$k <- 1L
b <- both(cl, 1, cluster_col = "k")
rec("FX-16", "RA-19 stride diverso per asse (x 2, y 3)", desc(b$o), desc(b$n),
    if (inherits(b$n$v, "err")) "errore esplicito invariato" else "NON atteso")
mm <- matrix(1L, 6, 6); mm[3, 3] <- 0L; mm[4, 4] <- 0L; i <- which(mm != 0, arr.ind = TRUE)
cl <- data.frame(x = i[, 2], y = i[, 1], k = mm[i])
b <- both(cl, 1, min_region_area_um2 = 0, simplify_tol_um = 0, cluster_col = "k")
rec("FX-17", "RA-19 due buchi diagonali, esatto (tol 0)", desc(b$o), desc(b$n),
    if (same_df(b$o$v$region_df, b$n$v$region_df) && !length(b$n$w) && b$n$v$region_df$n_holes == 2) "invariato, nessun warning" else "NON atteso")
b <- both(cl, 1, min_region_area_um2 = 0, cluster_col = "k")
rec("FX-18", "RA-19/RA-20 due buchi diagonali a 1 µm, tol 0.5 (= mezzo pixel)", desc(b$o), desc(b$n),
    if (same_df(b$o$v$region_df, b$n$v$region_df) && has_w(b$n, "simplify")) "invariato + warning di tolleranza (soglia >=, caso limite)" else "NON atteso",
    "tol = e/2 esatto: il warning scatta per '>='; a 1 µm la semplificazione sposta davvero la geometria (RA-21: I2, adv_pinch_hole)")
# ---- RA-20: semplificazione con pixel piccolo ------------------------------------------------------
d <- expand.grid(x = 1:121, y = 1:121); d <- d[(d$x - 61)^2 + (d$y - 61)^2 <= 60^2, ]; d$k <- 1L
b <- both(d, 0.25, cluster_col = "k")
rec("FX-19", "RA-20 disco r 60 px a 0.25 µm/px, tol 0.5", desc(b$o), desc(b$n),
    if (has_w(b$n, "simplify") && same_df(b$o$v$region_df, b$n$v$region_df)) "warning atteso presente, valori invariati" else "NON atteso")
ln <- data.frame(x = 1:800, y = 1L, k = 1L)
b <- both(ln, 0.25, min_region_area_um2 = 0, cluster_col = "k")
rec("FX-20", "RA-20 linea 800x1 px a 0.25 µm/px, tol 0.5", desc(b$o), desc(b$n),
    if (has_w(b$n, "simplify")) "warning atteso presente" else "NON atteso", "L'area resta 25 µm² invece di 50 (warning, non correzione)")
b <- both(ln, 0.25, min_region_area_um2 = 0, cluster_col = "k", simplify_tol_um = 0.12)
rec("FX-21", "RA-20 linea a 0.25 µm/px con tol 0.12 (< mezzo pixel)", desc(b$o), desc(b$n),
    if (!length(b$n$w) && abs(area_tot(b$n) - 50) < 1e-9) "nessun warning, area esatta" else "NON atteso")
b <- both(d, 1, cluster_col = "k", simplify_tol_um = 0.49)
rec("FX-22", "RA-20 disco a 1 µm/px con tol 0.49 (appena sotto mezzo pixel)", desc(b$o), desc(b$n),
    if (!length(b$n$w)) "nessun warning" else "NON atteso",
    sprintf("errore di area %.3g %% (diff. simmetrica non misurata qui)", 100 * (area_tot(b$n) / (nrow(d)) - 1)))
TT <- do.call(rbind, ROWS)
saveRDS(TT, file.path(REV, "10_fixcheck_adv.rds"))

# ---- regressione su input reali/sintetici/avversari del test --------------------------------------
idx <- read.csv("/home/user/2025.geo_spatialtrans/results/S1.1/S1.1_inputs.csv", stringsAsFactors = FALSE)
ids <- idx$id[idx$group != "perf"]
reg <- mclapply(ids, function(id) {
  o <- readRDS(file.path(IN, paste0(id, ".rds")))
  out <- list()
  for (cfg in list(c(min = 100, tol = 0.5), c(min = 0, tol = 0))) {
    a <- run(OLD, o$clust, o$pixel_size_um, min_region_area_um2 = cfg[["min"]], cluster_col = o$cluster_col, simplify_tol_um = cfg[["tol"]])
    n <- run(NEW, o$clust, o$pixel_size_um, min_region_area_um2 = cfg[["min"]], cluster_col = o$cluster_col, simplify_tol_um = cfg[["tol"]])
    ia <- a$v$info; inn <- n$v$info; ia$elapsed_s <- inn$elapsed_s <- NULL; ty_o <- typeof(ia$cluster_ids); ty_n <- typeof(inn$cluster_ids); if (is.numeric(ia$cluster_ids)) { ia$cluster_ids <- as.numeric(ia$cluster_ids); inn$cluster_ids <- as.numeric(inn$cluster_ids) }; src <- inn$stride_source; inn$stride_source <- NULL
    out[[length(out) + 1L]] <- data.frame(id = id, group = o$group, px = o$pixel_size_um, cfg = sprintf("min%g_tol%g", cfg[["min"]], cfg[["tol"]]),
      identical_region_df = same_df(a$v$region_df, n$v$region_df), identical_excluded = same_df(a$v$excluded_df, n$v$excluded_df),
      identical_polys = identical(st_as_binary(a$v$region_polygons), st_as_binary(n$v$region_polygons)),
      identical_info = identical(ia, inn), cluster_id_type_old = ty_o, cluster_id_type_new = ty_n, stride = inn$stride, stride_source = src,
      warn_old = length(a$w), warn_new = paste(sub("^extract_regions\\(\\): ", "", n$w), collapse = " || "))
  }
  do.call(rbind, out)
}, mc.cores = 12, mc.preschedule = FALSE)
bad <- vapply(reg, function(z) !is.data.frame(z), TRUE); if (any(bad)) print(reg[bad])
REG <- do.call(rbind, reg[!bad])
write.csv(REG, file.path(REV, "10_fixcheck_regression.csv"), row.names = FALSE)

# ---- dove scatterebbe il warning di stride su sottoinsiemi realistici (un cluster alla volta) -----
g <- function(v) { u <- sort(unique(v)); if (length(u) < 2) return(NA_integer_); d <- unique(diff(u)); r <- 0; for (a in d) { b <- r; while (a > 0) { t <- b %% a; b <- a; a <- t }; r <- b }; r }
sub <- do.call(rbind, lapply(ids[grepl("^real_|^roi_", ids)], function(id) {
  o <- readRDS(file.path(IN, paste0(id, ".rds"))); cl <- o$clust; k <- cl[[o$cluster_col]]
  do.call(rbind, lapply(split(seq_len(nrow(cl)), k), function(ii) {
    sx <- g(cl$x[ii]); sy <- g(cl$y[ii]); s <- if (is.na(sx) && is.na(sy)) 1 else if (is.na(sx)) sy else if (is.na(sy)) sx else if (sx != sy) -1 else sx
    data.frame(id = id, n = length(ii), s = s)
  }))
}))
write.csv(sub, file.path(REV, "10_fixcheck_cluster_subsets.csv"), row.names = FALSE)
cat(sprintf("SUBSETS: %d sottoinsiemi per cluster; stride stimato >1: %d; errore assi diversi: %d\n", nrow(sub), sum(sub$s > 1), sum(sub$s == -1)))

