# R/testing/test_S2.8.R - S2.8: batteria del consolidamento del codice dello Step 2 (impronta della ripresa con codice e profili, nomi della
# sovra-generazione, versione; BL-247). Cancello G-1 di docs/steps/S2.8.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.8.R [gruppi]
# Stampa PASS/FAIL per asserzione; esce con 1 se una fallisce. Variabile S28_RDIR: cartella dei sorgenti (mutanti).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
RDIR <- Sys.getenv("S28_RDIR", "R")
for (f in sort(list.files(RDIR, pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
GR <- commandArgs(TRUE); GR <- if (length(GR)) strsplit(GR[1], ",")[[1]] else c("a", "b", "c", "d")
nf <- 0L; np <- 0L
ok <- function(name, cond) { cond <- isTRUE(cond); cat(sprintf("%s  %s\n", if (cond) "PASS" else "FAIL", name)); if (cond) np <<- np + 1L else nf <<- nf + 1L }
err_of <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))
mk_layer <- function(preset, side, seed) {
  idx <- which(matrix(TRUE, side, side), arr.ind = TRUE); cl <- data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L)
  suppressWarnings(place_cells(cl, preset, pixel_size_um = 2, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                    small_components = "exclude", random_seed = seed, types = "multi", size_arm = "recentred", spatial = "uniform",
                    size_resid = "gaussian", strict = FALSE, verbose = FALSE))
}
L4 <- mk_layer("lymphoid", 120L, 104L)
PR <- load_cell_profiles("A4"); DF <- visium_hd_depth_factor("A4")
run <- function(dir, profiles = PR, ...) simulate_visium_hd_pieces(L4, seed = 7L, piece_um = 112, frame_open = TRUE, profiles = profiles,
                                                                    depth_factor = DF, out_dir = dir, ...)
# --- (a) impronta della ripresa
if ("a" %in% GR) {
  d1 <- tempfile("s28_a_"); o1 <- run(d1); fp <- readRDS(file.path(d1, "run.rds"))
  ok("a1 l'impronta contiene l'impronta del codice e quella dei profili",
     all(c("code", "profiles") %in% names(fp)) && identical(fp$code, .geo_code_hash("simulate_visium_hd_pieces")) &&
     identical(fp$profiles, .geo_object_hash(PR)) && grepl("^[0-9a-f]{32}$", fp$code) && grepl("^[0-9a-f]{32}$", fp$profiles))
  fs <- list.files(d1, pattern = "^(in|s0|A|B)_[0-9]+\\.rds$", full.names = TRUE); mt <- file.mtime(fs)
  Sys.sleep(1.1); o1b <- run(d1)
  ok(sprintf("a2 ripresa con codice e argomenti uguali: file dei pezzi riusati (%d) e uscite identiche", length(fs)),
     length(fs) > 0 && identical(file.mtime(fs), mt) && identical(o1$counts_2um, o1b$counts_2um) && identical(o1$truth$cells, o1b$truth$cells))
  orig <- .p27_tab
  .p27_tab <- function(key, w, n) { v <- numeric(n); if (length(key) > 0L) { s <- rowsum(as.numeric(w), key); v[as.integer(rownames(s))] <- s[, 1] }; v }
  e3 <- err_of(run(d1)); .p27_tab <- orig
  ok("a3 ripresa dopo una modifica di una funzione dello Step 2 (stesso risultato, codice diverso): si ferma e nomina il codice",
     grepl("esecuzione diversa", e3) && grepl("\\bcode\\b", e3) && !grepl("profiles", e3))
  zt <- .s23_ztab; .s23_ztab <- c(zt[1] - 1e-9, zt[-1]); e3b <- err_of(run(d1)); .s23_ztab <- zt
  ok("a3b ripresa dopo una modifica di una costante letta dal codice (tabella .s23_ztab): si ferma e nomina il codice",
     grepl("esecuzione diversa", e3b) && grepl("\\bcode\\b", e3b))
  PR2 <- PR; PR2$P[1, 1] <- PR2$P[1, 1] + 1
  e4 <- err_of(run(d1, profiles = PR2))
  ok("a4 ripresa con profili diversi (un conteggio in piu'): si ferma e nomina i profili, non il codice",
     grepl("esecuzione diversa", e4) && grepl("profiles", e4) && !grepl("\\bcode\\b", e4))
  fo <- fp; fo$code <- NULL; fo$profiles <- NULL; d5 <- tempfile("s28_a5_"); dir.create(d5); saveRDS(fo, file.path(d5, "run.rds"))
  e5 <- err_of(run(d5))
  ok("a5 cartella di un'esecuzione di S2.7 (impronta senza codice e profili): si ferma", grepl("esecuzione diversa", e5) && grepl("code", e5) && grepl("profiles", e5))
  h0 <- .geo_code_hash("simulate_visium_hd_pieces")
  txt <- paste(deparse(orig), collapse = "\n"); txt <- sub("\\{", "{ # solo un commento\n", txt)
  .p27_tab <- eval(parse(text = txt, keep.source = TRUE)); h1 <- .geo_code_hash("simulate_visium_hd_pieces"); .p27_tab <- orig
  ok("a6 un commento in piu' (sorgente conservata) non cambia l'impronta del codice", identical(h0, h1))
  pn <- sub("^## (<costante> )?", "", names(.geo_code_parts("simulate_visium_hd_pieces")))
  need <- c("simulate_visium_hd_pieces", ".p27_tab", ".s23_ztab", ".p27_stage0", ".p27_stageA", ".p27_stageB", ".p27_assemble", ".p27_place_tile", ".p27_edt_window",
            "simulate_stripes", "stripe_ceiling", "capture_params", "transcript_params", "load_cell_profiles", "hd_grid_geometry", ".s23_edt")
  ok(sprintf("a7 l'impronta copre le funzioni dei passi chiamate dalla funzione a pezzi (%d funzioni; mancanti: %s)", length(pn),
             paste(setdiff(need, pn), collapse = ", ")), all(need %in% pn))
  ok("a8 l'impronta del motore di S2.5 copre i quattro passi", all(c("simulate_visium_hd", "build_hd_grid", "simulate_cell_profiles",
       "place_and_move_transcripts", "capture_and_bin", "simulate_stripes") %in% sub("^## (<costante> )?", "", names(.geo_code_parts("simulate_visium_hd")))))
  L4b <- L4; ly_ <- L4b$cell_df$cell_type; sw <- which(ly_ %in% unique(ly_)[1:2]); L4b$cell_df$cell_type[sw] <- rev(ly_[sw])
  e9 <- err_of(simulate_visium_hd_pieces(L4b, seed = 7L, piece_um = 112, frame_open = TRUE, profiles = PR, depth_factor = DF, out_dir = d1))
  ok(sprintf("a9 ripresa con il livello cellulare cambiato (tipi di %d cellule scambiati, stessi centri): si ferma e nomina il livello", sum(ly_[sw] != rev(ly_[sw]))),
     sum(ly_[sw] != rev(ly_[sw])) > 0 && grepl("esecuzione diversa", e9) && grepl("\\blayer\\b", e9))
  gl <- lapply(get0_parts <- names(.geo_code_parts("simulate_visium_hd_pieces")), function(n) NULL)
  fnames <- sub("^## ", "", grep("^## [^<]", get0_parts, value = TRUE))
  G <- lapply(fnames, function(n) codetools::findGlobals(get(n), merge = FALSE)); af <- unique(unlist(lapply(G, `[[`, "functions"))); av <- unique(unlist(lapply(G, `[[`, "variables")))
  cand_call <- setdiff(intersect(c("nrow", "sum", "length", "seq_len", "min", "round", "floor", "paste0"), af), av)[1]
  cand_val <- intersect(c("c", "max", "rbind", "sum"), av)[1]
  assign(cand_call, 5, envir = globalenv()); h2 <- .geo_code_hash("simulate_visium_hd_pieces"); rm(list = cand_call, envir = globalenv())
  assign(cand_val, 5, envir = globalenv()); h3 <- .geo_code_hash("simulate_visium_hd_pieces"); rm(list = cand_val, envir = globalenv())
  ok(sprintf("a10 un oggetto dell'utente con il nome di una funzione di base usata solo in chiamata (%s) non cambia l'impronta; con il nome di una funzione passata come valore (%s, che il codice userebbe davvero) la cambia",
             cand_call, cand_val), !is.na(cand_call) && !is.na(cand_val) && identical(h0, h2) && !identical(h0, h3))
  sc <- tempfile(fileext = ".R")
  writeLines(c('.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")', 'suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })',
               sprintf('for (f in sort(list.files("%s", pattern = "^(04b|05b).*\\\\.R$", full.names = TRUE))) source(f)', RDIR),
               'cat(.geo_code_hash("simulate_visium_hd_pieces"), .geo_object_hash(load_cell_profiles("A4")), "\\n")'), sc)
  hh <- vapply(c("C.UTF-8", "C"), function(lc) trimws(paste(system2("Rscript", c("--vanilla", sc), stdout = TRUE, env = c(paste0("LC_ALL=", lc), "R_LIBS_SITE=/nonexistent", "R_LIBS_USER=/nonexistent")), collapse = "")), character(1))
  ok(sprintf("a11 impronte del codice e dei profili uguali in due processi nuovi con lingua C.UTF-8 e C e in questo processo (%s)", paste(substr(hh, 1, 8), collapse = " / ")),
     length(unique(hh)) == 1L && identical(strsplit(hh[[1]], " ")[[1]], c(h0, .geo_object_hash(PR))))
  unlink(c(d1, d5, sc), recursive = TRUE)
}
# --- (b) nomi della sovra-generazione
if ("b" %in% GR) {
  o <- simulate_visium_hd_pieces(L4, seed = 7L, piece_um = 112, frame_open = TRUE); st <- o$truth$layers$stripes
  Ms <- stripe_ceiling(capture_params("A4"))
  ok(sprintf("b1 truth$layers$stripes$M = metadata$overgeneration = tetto M* (%.4f)", Ms), identical(st$M, o$metadata$overgeneration) && identical(st$M, Ms))
  ok(sprintf("b2 M_raw = C = massimo grezzo dei fattori (%.4f) > M*: quadratini tagliati %d", st$M_raw, o$metadata$n_bins_over_ceiling),
     identical(st$M_raw, max(o$bins$c_raw)) && identical(st$C, st$M_raw) && st$M_raw > Ms && o$metadata$n_bins_over_ceiling == sum(o$bins$c_raw > Ms))
  ok("b3 bins$c = min(bins$c_raw, M*), esatto", identical(o$bins$c, pmin(o$bins$c_raw, Ms)))
  o1s <- simulate_visium_hd_pieces(L4, seed = 7L, piece_um = 112, frame_open = TRUE, layers = list(stripes = FALSE))
  rg <- sum(o$truth$cells$n_transcripts_generated) / sum(o1s$truth$cells$n_transcripts_generated)
  ok(sprintf("b4 profili generati con la profondita' per M*: trascritti generati con le strisce / senza = %.4f contro M* = %.4f (entro il 3 %%)", rg, Ms),
     abs(rg / Ms - 1) < 0.03 && isTRUE(all.equal(o$truth$cells$n_expected_captured, o$truth$cells$n_transcripts_generated / Ms)))
  lyr <- list(stripes = FALSE)
  oo <- simulate_visium_hd_pieces(L4, seed = 7L, piece_um = 112, frame_open = TRUE, layers = lyr)
  ok("b5 senza strisce: nessuno strato delle strisce, sovra-generazione 1, c = c_raw = 1", is.null(oo$truth$layers$stripes) && identical(oo$metadata$overgeneration, 1) &&
     all(oo$bins$c == 1) && all(oo$bins$c_raw == 1))
  s5 <- simulate_visium_hd(L4, seed = 7L, frame_open = TRUE)
  ok("b6 motore di S2.5: truth$layers$stripes$M = metadata$overgeneration = massimo di bins$c (stesso significato: sovra-generazione usata)",
     identical(s5$truth$layers$stripes$M, s5$metadata$overgeneration) && identical(s5$metadata$overgeneration, max(s5$bins$c)))
}
# --- (c) versione e impronta del codice nei metadati
if ("c" %in% GR) {
  v <- paste("geo.spatialtrans", unname(read.dcf("DESCRIPTION", fields = "Version")[1, 1]))
  o <- simulate_visium_hd_pieces(L4, seed = 7L, piece_um = Inf, frame_open = TRUE, layers = list(stripes = FALSE))
  s5 <- simulate_visium_hd(L4, seed = 7L, frame_open = TRUE, layers = list(stripes = FALSE))
  ok(sprintf("c1 metadata$package_version da DESCRIPTION nei due motori (%s)", v), identical(o$metadata$package_version, v) && identical(s5$metadata$package_version, v))
  ok("c2 impronta del codice nei metadati dei due motori", identical(o$metadata$code_hash, .geo_code_hash("simulate_visium_hd_pieces")) &&
     identical(s5$metadata$code_hash, .geo_code_hash("simulate_visium_hd")) && !identical(o$metadata$code_hash, s5$metadata$code_hash))
  ok("c3 impronta dei profili nei metadati della funzione a pezzi", identical(o$metadata$profiles_hash, .geo_object_hash(load_cell_profiles("A4"))))
}
# --- (d) test aggiunti dopo la mappa del documento di progetto (dichiarati nel file dello step)
if ("d" %in% GR) { if (file.exists("R/testing/test_S2.8_d.R")) source("R/testing/test_S2.8_d.R", local = TRUE) else ok("d il file del gruppo d esiste", FALSE) }
cat(sprintf("\n%d PASS, %d FAIL\n", np, nf))
quit(status = if (nf) 1L else 0L)
