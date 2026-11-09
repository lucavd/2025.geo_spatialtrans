# R/testing/test_S2.5.R - S2.5: batteria di test di simulate_visium_hd(), visium_hd_chunk_size(), visium_hd_archetype() (cancello G-1, parti a-f).
# Pre-registrazione: docs/steps/S2.5.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S25_RSRC=<cartella con R/>] Rscript --vanilla R/testing/test_S2.5.R
# I mutanti sono copie modificate del codice (tools/S2.5_mutants.sh) passate con S25_RSRC.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
SRC <- Sys.getenv("S25_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
err <- function(expr) inherits(tryCatch({ expr; NULL }, error = function(e) e), "error")

mk_layer <- function(nx, ny, a, b, hole = NULL, preset = "epithelial") {
  sq <- function(x0, y0, s) sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
  xy <- expand.grid(i = seq_len(nx) - 1L, j = seq_len(ny) - 1L)
  cx <- xy$i * a; cy <- xy$j * a
  if (!is.null(hole)) { k <- !(cx >= hole[1] & cx + a <= hole[3] & cy >= hole[2] & cy + a <= hole[4]); cx <- cx[k]; cy <- cy[k] }
  Tt <- sf::st_sfc(lapply(seq_along(cx), function(k) sq(cx[k], cy[k], a)))
  Nn <- sf::st_sfc(lapply(seq_along(cx), function(k) sq(cx[k] + (a - b) / 2, cy[k] + (a - b) / 2, b)))
  R <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(nx * a, 0), c(nx * a, ny * a), c(0, ny * a), c(0, 0)))))
  if (!is.null(hole)) R <- sf::st_difference(R, sf::st_sfc(sq(hole[1], hole[2], hole[3] - hole[1])))
  n <- length(cx)
  list(cell_df = data.frame(cell_id = 100L + seq_len(n), cell_type = rep(c("x", "y"), length.out = n), stringsAsFactors = FALSE),
       cell_territories = Tt, cell_nuclei = Nn, region_polygons = R,
       metadata = list(pixel_size_um = 1, random_seed = 1, tissue_preset_used = preset, package_version = "test"))
}
# profili costruiti: due tipi, tre sottotipi, 6 geni; totali per tipo da una distribuzione nota
mk_profiles <- function(N = 400, depth = 1) {
  P <- rbind(x1 = c(.4, .3, .1, .1, .05, .05), x2 = c(.05, .05, .1, .1, .3, .4), y1 = rep(1 / 6, 6))
  colnames(P) <- sprintf("G%d", 1:6)
  list(archetype = "A1", genes = data.frame(gene_id = colnames(P), gene_name = colnames(P)),
       subtypes = data.frame(subtype = rownames(P), cell_type = c("x", "x", "y"), share = c(.6, .4, 1)), P = P, phi = rep(0, 6),
       totals = data.frame(q = c(0, 1), x = c(N * 0.8, N * 1.2), y = c(N * 0.5, N * 1.5)), totals_subtype = NULL, noise = "M",
       depth_factor = depth, meta = list(totals_by = "type", winner = "test"))
}
L <- mk_layer(10, 8, 10, 4, hole = c(40, 30, 60, 50)); H <- build_hd_grid(L, margin_um = 16); PR <- mk_profiles()
ALL <- list(stripes = TRUE, background = TRUE, field = TRUE, second_tail = TRUE)
CS <- c(stripes = 11, profiles = 22, place = 33, capture = 44)
manual <- function(L, H, PR, A, cs, ly = ALL, chunk = 7L, frame_open = FALSE, depth = PR$depth_factor) {
  cp <- capture_params(A); p <- transcript_params(A); p$w2 <- if (ly$second_tail) cp$w2 else 0; p$ell2 <- cp$ell2
  St <- if (ly$stripes) simulate_stripes(H, cp, seed = cs[["stripes"]]) else NULL; M <- if (is.null(St)) 1 else St$M
  pr <- simulate_cell_profiles(L$cell_df[, c("cell_id", "cell_type")], PR, seed = cs[["profiles"]], depth_factor = depth * M)
  tx <- place_and_move_transcripts(L, H, pr, p, seed = cs[["place"]], field = ly$field, frame_open = frame_open, chunk_size = chunk)
  list(cap = capture_and_bin(tx, H, St, list(b = cp$b), seed = cs[["capture"]], overgeneration = M, background = ly$background), pr = pr, tx = tx, M = M)
}
same_out <- function(o, m) identical(o$counts_2um, m$cap$counts_2um) && identical(o$truth$origin, m$cap$truth$origin) &&
  identical(o$truth$background, m$cap$truth$background) && identical(o$truth$cell_counts, m$pr$counts)

# ================= (a) equivalenza con la catena a mano
for (A in c("A2", "A4")) {
  o <- simulate_visium_hd(L, archetype = A, profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H)
  m <- manual(L, H, PR, A, CS)
  ck(sprintf("a1 %s orchestratore = catena a mano (stessi semi e blocchi)", A), same_out(o, m), sprintf("M %.3f", o$metadata$overgeneration))
}
o2 <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H, frame_open = TRUE)
m2 <- manual(L, H, PR, "A2", CS, frame_open = TRUE)
ck("a2 frame_open passato alla posizione e allo spostamento", same_out(o2, m2) && !identical(o2$counts_2um, simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H)$counts_2um))
oH <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L)
ck("a3 senza hd: griglia costruita con build_hd_grid() (margine 280 um)", identical(oH$metadata$grid, build_hd_grid(L)$grid) && sum(oH$truth$cells$n_captured) > 0)

# ================= (b) identita' delle uscite
o <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H)
g <- H$grid; nb <- g$n_row * g$n_col; tr <- o$truth; cl <- tr$cells
ck("b1 counts_8um = somma a 8 um di counts_2um (S2.1)", identical(o$counts_8um, sum_to_8um(o$counts_2um, g)))
ck("b2 origin_8um = sum_to_8um(origin, by = cell_id)", { s <- sum_to_8um(tr$origin, g, cols = c("n", "n_born_nucleus"), by = "cell_id")
  isTRUE(all.equal(s, tr$origin_8um[, names(s)], check.attributes = FALSE)) && all(s$n == tr$origin_8um$n) })
ck("b3 overlap_8um = sum_to_8um(overlap, by = cell_id), frazioni su 64 um^2", { s <- sum_to_8um(tr$overlap, g, cols = c("area_territory_um2", "area_nucleus_um2"), by = "cell_id")
  isTRUE(all.equal(s$area_territory_um2, tr$overlap_8um$area_territory_um2, tolerance = 1e-12)) && identical(s$cell_id, tr$overlap_8um$cell_id) &&
  isTRUE(all.equal(tr$overlap_8um$frac_territory, tr$overlap_8um$area_territory_um2 / 64)) })
ck("b4 background_8um = sum_to_8um(background)", identical(as.numeric(sum_to_8um(tr$background, g, cols = "n")$n), as.numeric(tr$background_8um$n)))
ob <- numeric(nb); tb <- tapply(tr$origin$n, tr$origin$bin_id, sum); ob[as.integer(names(tb))] <- tb
bgv <- numeric(nb); bgv[tr$background$bin_id] <- tr$background$n
ck("b5 UMI per quadratino = origine + fondo", identical(as.numeric(Matrix::colSums(o$counts_2um)), ob + bgv), sprintf("fondo %d", sum(bgv)))
ck("b6 per cellula generati = catturati + scartati + persi; = colonne di cell_counts", all(cl$n_generated == cl$n_captured + cl$n_thinned + cl$n_lost) &&
     identical(as.numeric(cl$n_generated), as.numeric(Matrix::colSums(tr$cell_counts))) && identical(cl$cell_id, L$cell_df$cell_id) && identical(cl$cell_type, L$cell_df$cell_type))
ck("b7 arrivati per quadratino: somma = generati - persi; catturati <= arrivati", sum(tr$arrived_2um) == sum(cl$n_generated) - sum(cl$n_lost) && all(ob <= tr$arrived_2um))
ck("b8 atteso catturato = generati / M; fattore di profondita' per M nei profili", isTRUE(all.equal(cl$n_expected_captured, cl$n_transcripts_generated / o$metadata$overgeneration)) &&
     isTRUE(all.equal(tr$layers$depth$overgeneration, o$metadata$overgeneration)) && o$metadata$overgeneration > 1)
ck("b9 metadati: scala, semi, strati, blocchi", identical(o$metadata$pixel_size_um, 1) && identical(o$metadata$component_seeds, CS) && identical(o$metadata$chunk_size, 7L) &&
     isTRUE(all(unlist(o$metadata$layers))))

# ================= (c) strati spenti
oS <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H, layers = list(stripes = FALSE))
ck("c1 strisce spente: M = 1, c = 1, nessuno scartato, = catena a mano", oS$metadata$overgeneration == 1 && all(oS$bins$c == 1) && sum(oS$truth$cells$n_thinned) == 0 &&
     same_out(oS, manual(L, H, PR, "A2", CS, ly = modifyList(ALL, list(stripes = FALSE)))))
oB <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H, layers = list(background = FALSE))
ck("c2 fondo spento: nessun fondo; resto = uscita completa", nrow(oB$truth$background) == 0 && identical(oB$truth$origin, o$truth$origin))
oT <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H, layers = list(second_tail = FALSE))
ck("c3 seconda coda spenta = catena a mano con w2 = 0", oT$truth$layers$diffusion$w2 == 0 && same_out(oT, manual(L, H, PR, "A2", CS, ly = modifyList(ALL, list(second_tail = FALSE)))))
oF <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H, layers = list(field = FALSE))
ck("c4 campo spento = catena a mano senza campo", !oF$truth$layers$field$on && same_out(oF, manual(L, H, PR, "A2", CS, ly = modifyList(ALL, list(field = FALSE)))))
oD <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, chunk_size = 7L, hd = H, depth_factor = 2.5)
ck("c5 depth_factor sostituisce quello dei profili (per M)", same_out(oD, manual(L, H, PR, "A2", CS, depth = 2.5)) && oD$truth$layers$depth$depth_factor == 2.5)

# ================= (d) blocchi dai trascritti
ck("d1 visium_hd_chunk_size: 4e6 / 30000 = 133; tetto 2000; minimo 1; senza trascritti 2000",
   visium_hd_chunk_size(rep(30000, 10)) == 133L && visium_hd_chunk_size(rep(100, 10)) == 2000L && visium_hd_chunk_size(rep(1e8, 3)) == 1L && visium_hd_chunk_size(numeric(0)) == 2000L &&
   visium_hd_chunk_size(rep(1000, 5), chunk_transcripts = 5e4) == 50L)
ck("d2 chunk_transcripts non valido: errore", err(visium_hd_chunk_size(1, chunk_transcripts = 0)) && err(visium_hd_chunk_size(1, chunk_transcripts = NA)))
oA <- simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, hd = H, chunk_transcripts = 5000)
ck("d3 blocchi automatici = regola sui trascritti generati (sovra-generazione compresa)", oA$metadata$chunk_size == visium_hd_chunk_size(oA$truth$cells$n_transcripts_generated, 5000) &&
     same_out(oA, manual(L, H, PR, "A2", CS, chunk = oA$metadata$chunk_size)), sprintf("blocchi %d", oA$metadata$chunk_size))
ck("d4 chunk_size esplicito rispettato", simulate_visium_hd(L, archetype = "A2", profiles = PR, component_seeds = CS, hd = H, chunk_size = 3L)$metadata$chunk_size == 3L)

# ================= (e) semi, processi, generatore dell'utente
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(99); s0 <- .Random.seed
e1 <- simulate_visium_hd(L, seed = 2026, archetype = "A2", profiles = PR, hd = H, chunk_size = 7L)
ck("e1 generatore dell'utente intatto (tipo e stato)", identical(.Random.seed, s0) && identical(RNGkind()[1], "Mersenne-Twister"))
RNGkind("L'Ecuyer-CMRG"); set.seed(5); s1 <- .Random.seed
e1b <- simulate_visium_hd(L, seed = 2026, archetype = "A2", profiles = PR, hd = H, chunk_size = 7L)
ck("e2 stesso seme stesso risultato (anche con un altro generatore dell'utente, ripristinato)", identical(e1$counts_2um, e1b$counts_2um) && identical(.Random.seed, s1) && identical(RNGkind()[1], "L'Ecuyer-CMRG"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
exp_seeds <- local({ set.seed(2026); s <- sample.int(1e8, 4L); c(stripes = s[1], profiles = s[2], place = s[3], capture = s[4]) })
ck("e3 semi dei passi = 4 estrazioni Mersenne-Twister dal seme dell'utente", identical(e1$metadata$component_seeds, exp_seeds))
e2 <- simulate_visium_hd(L, seed = 2027, archetype = "A2", profiles = PR, hd = H, chunk_size = 7L)
ck("e4 seme diverso, risultato diverso; semi dei passi non consecutivi", !identical(e1$counts_2um, e2$counts_2um) && all(abs(diff(e1$metadata$component_seeds)) > 4))
e3 <- simulate_visium_hd(L, seed = 2026, archetype = "A2", profiles = PR, hd = H, chunk_size = 7L, workers = 2L)
ck("e5 1 e 2 processi identici", identical(e1$counts_2um, e3$counts_2um) && identical(e1$truth$origin, e3$truth$origin) && identical(e1$truth$background, e3$truth$background))

# ================= (f) errori espliciti e archetipo
Lp <- L; Lp$metadata$pixel_size_um <- NULL
ck("f1 scala mancante: errore", err(simulate_visium_hd(Lp, seed = 1, archetype = "A2", profiles = PR)))
ck("f2 seme mancante: errore", err(simulate_visium_hd(L, archetype = "A2", profiles = PR, hd = H)))
ck("f3 archetipo sconosciuto o non ricavabile: errore", err(simulate_visium_hd(L, seed = 1, archetype = "A9", profiles = PR, hd = H)) &&
     err(visium_hd_archetype(mk_layer(2, 2, 10, 4, preset = "boh"))))
ck("f4 archetipo dal preset", visium_hd_archetype(L) == "A1" && visium_hd_archetype(mk_layer(2, 2, 10, 4, preset = "columnar_muscle")) == "A6")
ck("f5 strati non validi: errore", err(simulate_visium_hd(L, seed = 1, archetype = "A2", profiles = PR, hd = H, layers = list(stripe = FALSE))) &&
     err(simulate_visium_hd(L, seed = 1, archetype = "A2", profiles = PR, hd = H, layers = list(stripes = NA))))
L2 <- mk_layer(12, 8, 10, 4)
ck("f6 griglia di un altro livello: errore", err(simulate_visium_hd(L2, seed = 1, archetype = "A2", profiles = PR, hd = H)))
ck("f7 chunk_size non valido: errore", err(simulate_visium_hd(L, seed = 1, archetype = "A2", profiles = PR, hd = H, chunk_size = 0)))
ck("f8 component_seeds incompleti: errore", err(simulate_visium_hd(L, archetype = "A2", profiles = PR, hd = H, component_seeds = c(stripes = 1, place = 2))))

cat(sprintf("\nTOTALE: %d PASS, %d FAIL\n", sum(unlist(res)), sum(!unlist(res))))
quit(status = if (all(unlist(res))) 0L else 1L)
