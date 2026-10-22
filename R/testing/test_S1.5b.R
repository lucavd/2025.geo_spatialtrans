# R/testing/test_S1.5b.R - S1.5b: test dei preset per tipo (cancello G-1; mutanti G-2). Pre-registrazione: docs/steps/S1.5b.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S15B_MUTANT=none|M1..M7] Rscript --vanilla R/testing/test_S1.5b.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "t_err"))
isok <- function(o) !inherits(o, "t_err"); msg <- function(o) if (isok(o)) "" else o$err
ARMS <- c("group", "recentred")

# ---- riferimenti letti dal file: righe B-186...B-208 ----
bio <- readLines("docs/BIO_REFERENCES.md"); bio <- bio[grepl("^\\| B-[0-9]{3} \\|", bio)]
cells <- function(l) trimws(strsplit(l, "|", fixed = TRUE)[[1]][-1])
B <- do.call(rbind, lapply(bio, function(l) { x <- cells(l); data.frame(id = x[1], arch = x[2], metric = x[3], value = x[4], dist = x[5], stringsAsFactors = FALSE) }))
B <- B[grepl("preset per tipo \\(S1.5b\\)", B$metric), ]
kv <- function(s, k) { m <- regmatches(s, regexpr(paste0("(^|; )", k, "=[^;]+"), s)); if (!length(m)) NA_real_ else as.numeric(sub(".*=", "", m)) }
B$type <- sub(".*\\*\\*(.+)\\*\\*.*", "\\1", B$metric)
B$fraction <- as.numeric(B$value)
for (k in c("gruppo", "rho", "beta", "sigma", "m_t", "m_sim", "frazione_follicolo", "frazione_paracorticale")) B[[k]] <- vapply(B$dist, kv, 0, k = k)
PRE <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
B$preset <- PRE[B$arch]
ck("T0 righe B-186...B-208 presenti, 23 tipi, 6 tessuti", nrow(B) == 23 && identical(B$id, sprintf("B-%d", 186:208)) && length(unique(B$arch)) == 6, nrow(B))
TT <- .TP_TABLE

# ---- T1: catalogo per tipo ricavato in modo indipendente dalle righe di BIO_REFERENCES ----
expect_cat <- function(p, arm, comp = "tissue") {
  b <- B[B$preset == p, ]; t <- TT[TT$preset == p, ]
  f <- switch(comp, tissue = b$fraction, follicle = b$frazione_follicolo, paracortex = b$frazione_paracorticale); f <- f / sum(f)
  aT <- sum(f * b$m_t)                                   # territorio medio della miscela (pesato per numero)
  dens <- t$density * aT / b$m_t                         # 1 / sum(f / dens) = D
  A0 <- 1e6 / t$density_model; At <- 1e6 / dens
  # rho_cat tale che log(rho^2 At) + beta (log A - log At) = log(rho_g^2 A0) + beta (log A - log A0) [+ beta log(m_t / m_sim)]
  lr2 <- log(b$rho^2) + (1 - b$beta) * (log(A0) - log(At))
  if (arm == "recentred") lr2 <- lr2 + b$beta * log(b$m_t / b$m_sim)
  data.frame(cell_type = b$type, density = dens, rho = sqrt(exp(lr2)), beta = b$beta, sigma = b$sigma, f = f, group = b$gruppo,
             min_dist = t$min_dist_um, phi = t$phi, pale = 1 - t$density_evident / t$density, stringsAsFactors = FALSE)
}
bad <- character(0)
for (p in tissue_presets()) for (arm in ARMS) {
  sp <- safe(tissue_preset_spec(p, types = "multi", size_arm = arm))
  if (!isok(sp)) { bad <- c(bad, paste(p, arm, msg(sp))); next }
  e <- expect_cat(p, arm); ct <- sp$cell_types[match(e$cell_type, sp$cell_types$cell_type), ]
  d <- c(density = max(abs(ct$density / e$density - 1)), rho = max(abs(ct$nucleus_to_eq_ratio / e$rho - 1)),
         beta = max(abs(ct$nucleus_area_slope - e$beta)), sigma = max(abs(ct$nucleus_area_sdlog - e$sigma)),
         f = max(abs(sp$composition$fraction[match(e$cell_type, sp$composition$cell_type)] - e$f)),
         group = max(abs(ct$size_group - e$group)), dmin = max(abs(ct$min_dist_um - e$min_dist)), phi = max(abs(ct$nucleus_nc_max - e$phi)),
         pale = max(abs(ct$pale_fraction - e$pale)))
  if (anyNA(d) || any(d > 1e-9) || nrow(sp$cell_types) != nrow(e)) bad <- c(bad, paste(p, arm, names(d)[is.na(d) | d > 1e-9]))
}
ck("T1 catalogo per tipo (densita', rho, beta, sigma, frazione, gruppo, distanza minima, phi, pallidi) uguale a quello ricavato da BIO_REFERENCES, 6 preset x 2 bracci", !length(bad), bad)

# ---- T2, T3: densita' della miscela e composizioni ----
bad2 <- character(0); bad3 <- character(0)
for (p in tissue_presets()) for (comp in c("tissue", if (p == "lymphoid") c("follicle", "paracortex"))) {
  sp <- safe(tissue_preset_spec(p, types = "multi", size_arm = "group", composition = comp))
  if (!isok(sp)) { bad2 <- c(bad2, paste(p, comp, msg(sp))); bad3 <- c(bad3, paste(p, comp)); next }
  mix <- 1 / sum(sp$composition$fraction / sp$cell_types$density[match(sp$composition$cell_type, sp$cell_types$cell_type)])
  D <- TT$density[TT$preset == p]
  if (abs(mix / D - 1) > 1e-9) bad2 <- c(bad2, sprintf("%s %s %.6f vs %d", p, comp, mix, D))
  b <- B[B$preset == p, ]; raw <- switch(comp, tissue = b$fraction, follicle = b$frazione_follicolo, paracortex = b$frazione_paracorticale)
  if (abs(sum(sp$composition$fraction) - 1) > 1e-12 || abs(sum(raw) - 1) > 1e-5) bad3 <- c(bad3, paste(p, comp, sum(raw)))
}
ck("T2 media armonica della miscela = densita' manuale del tessuto (B-043...B-048), anche nelle varianti del linfonodo", !length(bad2), bad2)
ck("T3 composizioni che sommano a 1 (righe di BIO_REFERENCES entro 1e-5, catalogo entro 1e-12)", !length(bad3), bad3)

# ---- T4: gruppi di taglia come la regola applicata (results/S1.5b) ----
cal <- read.csv("results/S1.5b/S1.5b_type_calibration.csv", stringsAsFactors = FALSE)
gB <- B$gruppo[match(paste(cal$archetype, cal$cell_type), paste(B$arch, B$type))]
ck("T4 gruppi di taglia di BIO_REFERENCES uguali a results/S1.5b/S1.5b_type_calibration.csv", nrow(cal) == 23 && all(gB == cal$size_group))
num_ok <- all(abs(B$rho[match(paste(cal$archetype, cal$cell_type), paste(B$arch, B$type))] / cal$rho_group - 1) < 5e-6) &&
  all(abs(B$fraction[match(paste(cal$archetype, cal$cell_type), paste(B$arch, B$type))] / cal$fraction_tissue - 1) < 5e-6)
ck("T4b rho dei gruppi e frazioni di BIO_REFERENCES uguali alla calibrazione (6 cifre)", num_ok)
ng <- tapply(cal$size_group, cal$archetype, function(g) length(unique(g)))
ck("T4c numero di gruppi per tessuto: 2, 2, 2, 1, 2, 1 (regola dei gruppi, docs/steps/S1.5b.md)", identical(as.integer(ng), c(2L, 2L, 2L, 1L, 2L, 1L)), paste(ng))

# ---- T5: il preset a un tipo non cambia ----
same <- all(vapply(tissue_presets(), function(p) identical(tissue_preset_spec(p)$cell_types, tissue_preset_spec(p, types = "single")$cell_types) &&
                     nrow(tissue_preset_spec(p)$cell_types) == 1L && is.null(tissue_preset_spec(p)$composition), TRUE))
ck("T5 types = \"single\" (predefinito durante il banco) restituisce il catalogo di un tipo di S1.5", same)

# ---- T6: place_cells con i preset per tipo ----
toy <- function(nx = 260, ny = 200) { cl <- expand.grid(x = 1:nx, y = 1:ny); cl$intensity_cluster <- ifelse(cl$x < nx * 0.45, 1L, 2L); cl }
cl0 <- toy(); bad6 <- character(0); bad6b <- character(0); bad6c <- character(0)
for (p in tissue_presets()) for (arm in ARMS) {
  L <- safe(suppressWarnings(place_cells(cl0, p, pixel_size_um = if (p %in% c("brain_cortex", "columnar_muscle", "stromal_rich")) 2 else 1,
                                          types = "multi", size_arm = arm, random_seed = 5, verbose = FALSE)))
  if (!isok(L)) { bad6 <- c(bad6, paste(p, arm, msg(L))); next }
  if (!L$metadata$validation$ok) bad6 <- c(bad6, paste(p, arm, "ERROR"))
  if (sum(L$region_df$n_failed) > 0) bad6b <- c(bad6b, paste(p, arm, "non piazzate", sum(L$region_df$n_failed)))
  comp <- L$region_composition
  for (r in seq_len(nrow(L$region_df))) {
    rr <- L$region_df[r, ]; n <- rr$n_cells; cdr <- L$cell_df[L$cell_df$region_id == rr$region_id, ]
    cr <- comp[as.character(comp$cluster_id) == as.character(rr$cluster_id), ]
    obs <- as.numeric(table(factor(as.character(cdr$cell_type), levels = cr$cell_type)))
    if (any(abs(obs - n * cr$fraction) >= 1 + 1e-9)) bad6c <- c(bad6c, paste(p, arm, rr$region_id))
  }
  md <- L$metadata
  if (!identical(md$preset_types, "multi") || !identical(md$size_arm, arm) || is.null(md$type_sources)) bad6 <- c(bad6, paste(p, arm, "metadati"))
}
ck("T6 place_cells sui 6 preset per tipo e 2 bracci: nessun errore, controlli senza ERROR, metadati del preset", !length(bad6), bad6)
ck("T6b nessuna cellula non piazzata (BL-055: distanza minima unica per tessuto)", !length(bad6b), bad6b)
ck("T6c conteggi per tipo e regione entro 1 dalla ripartizione n * frazione (arrotondamento di Madow)", !length(bad6c), bad6c)

# ---- T7, T8, T9: epitelio, varianti, errori ----
lab_ok <- TRUE
for (p in tissue_presets()) { ct <- tissue_preset_spec(p, types = "multi", size_arm = "group")$cell_types
  want <- ifelse(ct$cell_type == "epiteliale" & p %in% c("tumor_microenv", "stromal_rich"), "tumorale", NA_character_)
  lab_ok <- lab_ok && identical(unname(ct$epithelium_label), want) }
ck("T7 epitelio del tumore e dello stroma qualificato «tumorale» per istologia (DEC-014), nessuna qualifica altrove", lab_ok)
e1 <- safe(tissue_preset_spec("epithelial", types = "multi", size_arm = "group", composition = "follicle"))
f1 <- tissue_preset_spec("lymphoid", types = "multi", size_arm = "group", composition = "follicle")$composition
ck("T8 varianti del linfonodo: follicolo con B > T come le righe di BIO_REFERENCES; variante assente negli altri preset -> errore",
   !isok(e1) && f1$fraction[f1$cell_type == "linfociti B"] > f1$fraction[f1$cell_type == "linfociti T"])
ck("T9 types = \"multi\" senza size_arm o con size_conversion = \"none\" -> errore",
   !isok(safe(tissue_preset_spec("epithelial", types = "multi"))) && !isok(safe(tissue_preset_spec("epithelial", "none", types = "multi", size_arm = "group"))))

# ---- T10: due taglie nel cervello e relazione fra i bracci ----
g <- tissue_preset_spec("brain_cortex", types = "multi", size_arm = "group")$cell_types; rc <- tissue_preset_spec("brain_cortex", types = "multi", size_arm = "recentred")$cell_types
medA <- function(ct, t, A = 800) { k <- ct$cell_type == t; Ar <- 1e6 / ct$density[k]; ct$nucleus_to_eq_ratio[k]^2 * Ar * (A / Ar)^ct$nucleus_area_slope[k] }
ck("T10 cervello: a parita' di territorio (800 um2) il nucleo mediano dei neuroni e' oltre 2.1 volte quello degli oligodendrociti, in entrambi i bracci",
   medA(g, "neuroni") / medA(g, "oligodendrociti") > 2.1 && medA(rc, "neuroni") / medA(rc, "oligodendrociti") > 2.1,
   sprintf("%.2f %.2f", medA(g, "neuroni") / medA(g, "oligodendrociti"), medA(rc, "neuroni") / medA(rc, "oligodendrociti")))

cat(sprintf("RISULTATO test_S1.5b: %d/%d PASS\n", sum(unlist(res)), length(res)))
