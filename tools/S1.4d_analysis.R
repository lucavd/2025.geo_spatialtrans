# tools/S1.4d_analysis.R — S1.4d: cancelli G-3, G-4 e previsioni P-1...P-5, controprova CP-1, con le regole pre-registrate
# (docs/steps/S1.4d.md, commit 5118bce / eae40f6). Uso (root del repo): Rscript --vanilla tools/S1.4d_analysis.R
# Uscite: results/S1.4d/S1.4d_*.csv
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
OUT <- "/mnt/micron/geo_spatialtrans/S1.4d"; RES <- "results/S1.4d"; dir.create(RES, showWarnings = FALSE)
TIS <- c(A1 = "intestino", A2 = "tumore del colon", A3 = "stroma", A4 = "linfonodo", A5 = "cervello", A6 = "cuore")
rd <- function(d) { z <- lapply(sort(list.files(file.path(OUT, d), pattern = "\\.rds$", full.names = TRUE)), readRDS)
  nm <- unique(unlist(lapply(z, names)))                       # unione delle colonne (le unita' "fold" hanno le colonne eq_*)
  do.call(rbind, lapply(z, function(x) { for (m in setdiff(nm, names(x))) x[[m]] <- NA; x[, nm] })) }
# stima +- 2 SE; per le frazioni gli estremi si tagliano a [0, 1] (definizione operativa scritta prima del calcolo)
judge <- function(est, se, L, U, frac = TRUE) {
  lo <- est - 2 * se; hi <- est + 2 * se; if (frac) { lo <- max(0, lo); hi <- min(1, hi) }
  eps <- 1e-12
  v <- if (lo >= L - eps && hi <= U + eps) "confermata" else if (hi < L - eps || lo > U + eps) "smentita" else "indeterminata"
  data.frame(est = est, se = se, lo2 = lo, hi2 = hi, int_lo = L, int_hi = U, verdict = v)
}
gates <- list(); G <- function(id, check, value, thr, pass) gates[[length(gates) + 1]] <<- data.frame(gate = id, check = check, value = as.character(value), threshold = thr, pass = pass)

# ---- G-3 invarianti -------------------------------------------------------------------------------------------------
gen <- rd("gen"); nul <- rd("null"); pf <- readRDS(file.path(OUT, "perf", "perf_S.rds")); rep <- rd("repro")
inv_cols <- c("n", "c3_geos", "c3_sf", "no_pairs", "no_max", "area_rel_max_uncapped", "area_rel_max_capped", "n_invalid", "n_empty", "n_multipart", "n_clipped", "n_capped", "max_iter")
sets <- list(primary_tissue = gen[gen$method != "stardist_he" & gen$phimode == "tissue", ], primary_fold = gen[gen$phimode == "fold", ],
             stardist = gen[gen$method == "stardist_he", ], null = nul, section_S = pf$inv)
g3 <- do.call(rbind, lapply(names(sets), function(s) { x <- sets[[s]]
  data.frame(set = s, units = nrow(x), cells = sum(x$n), c3_geos = sum(x$c3_geos), c3_sf = sum(x$c3_sf), overlap_pairs = sum(x$no_pairs),
             overlap_max = max(x$no_max), area_rel_max = max(c(x$area_rel_max_uncapped, x$area_rel_max_capped)), invalid = sum(x$n_invalid),
             empty = sum(x$n_empty), multipart = sum(x$n_multipart), clipped = sum(x$n_clipped), capped = sum(x$n_capped), max_iter = max(x$max_iter)) }))
write.csv(g3, file.path(RES, "S1.4d_g3_invariants.csv"), row.names = FALSE)
G("G-3", "C3 (geos), nuclei fuori dal territorio", sum(g3$c3_geos), "0", sum(g3$c3_geos) == 0)
G("G-3", "C3 (sf via WKB), nuclei fuori dal territorio", sum(g3$c3_sf), "0", sum(g3$c3_sf) == 0)
G("G-3", "coppie sovrapposte > 1e-9 µm²", sum(g3$overlap_pairs), "0", sum(g3$overlap_pairs) == 0)
G("G-3", "area relativa max |A - A*| / A*", signif(max(g3$area_rel_max), 3), "<= 1e-6", max(g3$area_rel_max) <= 1e-6)
G("G-3", "geometrie invalide / vuote", paste(sum(g3$invalid), sum(g3$empty)), "0 / 0", sum(g3$invalid) == 0 && sum(g3$empty) == 0)
G("G-3", "nulli rigenerati = S1.4 (n, seme)", sum(nul$same_as_S14), "60/60", all(nul$same_as_S14) && nrow(nul) == 60)
G("G-3", "riproducibilita' bit per bit (md5)", sum(rep$md5_gen == rep$md5_repro), "6/6", nrow(rep) == 6 && all(rep$md5_gen == rep$md5_repro))

# ---- G-4 equivalenza ------------------------------------------------------------------------------------------------
fo <- gen[gen$phimode == "fold", ]
g4a <- do.call(rbind, lapply(split(fo, fo$archetype), function(x) data.frame(archetype = x$archetype[1], units = nrow(x), cells = sum(x$n),
  dA_rel_max = max(x$eq_dA_rel_max), sym_rel_max = max(x$eq_sym_rel_max), cells_sym_gt1e6 = sum(x$eq_sym_rel_gt1e6),
  cells_beyond_grid_bound = sum(x$eq_beyond_bound), units_same_ncap = sum(x$eq_ncap_same))))
write.csv(g4a, file.path(RES, "S1.4d_g4a_equivalence.csv"), row.names = FALSE)
G("G-4a", "forma pre-registrata: celle con differenza simmetrica / A* > 1e-6", sum(g4a$cells_sym_gt1e6), "0", sum(g4a$cells_sym_gt1e6) == 0)
G("G-4a", "stesse celle al tetto del prototipo (unita' x seme)", sum(g4a$units_same_ncap), "300/300", sum(g4a$units_same_ncap) == 300 && nrow(fo) == 300)
G("G-4a*", "spiegazione: |dA| / A* max (accordo d'area)", signif(max(g4a$dA_rel_max), 3), "<= 1e-6", max(g4a$dA_rel_max) <= 1e-6)
G("G-4a*", "spiegazione: celle oltre il limite del rumore di griglia", sum(g4a$cells_beyond_grid_bound), "0", sum(g4a$cells_beyond_grid_bound) == 0)
lg <- rd("legacy"); write.csv(lg, file.path(RES, "S1.4d_g4b_legacy.csv"), row.names = FALSE)
G("G-4b", "modo scale = pacchetto db4afea (cell_df e WKB)", sum(lg$identical_db4afea), "330/330", nrow(lg) == 330 && all(lg$identical_db4afea))
G("G-4b", "modo scale = righe V-c di S1.4b (C-0)", sum(lg$c0_identical_S14b), "330/330", nrow(lg) == 330 && all(lg$c0_identical_S14b))

# ---- previsioni -----------------------------------------------------------------------------------------------------
ci <- read.csv("results/R2c/R2c_posthoc_ci.csv"); ci <- ci[ci$ci == "due stadi", ]; rownames(ci) <- ci$archetype
w <- read.csv(file.path(RES, "S1.4d_windows.csv"))
k2r <- read.csv("results/R2c/R2c_k2_windows.csv"); pr <- w[w$source == "K2_proto", ]
m <- merge(k2r, pr, by = c("win_id", "seed"), suffixes = c(".r2c", ".d"))
G("controllo", "codice delle finestre = R2c (prototipo K-2: n e contatti identici)", sum(m$n.r2c == m$n.d & m$n_contact.r2c == m$n_contact.d), "240/240",
  nrow(m) == 240 && all(m$n.r2c == m$n.d & m$n_contact.r2c == m$n_contact.d))
s1 <- function(src, A) { x <- w[w$source == src & w$archetype == A, ]; ps <- tapply(x$s1, x$seed, mean)   # finestre a peso uguale per seme
  cbind(judge(mean(ps), sd(ps) / sqrt(length(ps)), ci[A, "lo"], ci[A, "hi"]), n_units = length(ps), n_windows = length(unique(x$win_id))) }
P <- list()
for (A in names(TIS)) P[[length(P) + 1]] <- cbind(id = "P-1", archetype = A, tissue = TIS[[A]], s1("K2_tissue", A))
P[[length(P) + 1]] <- cbind(id = "P-4", archetype = "A4", tissue = TIS[["A4"]], s1("K2_stardist", "A4"))
d_a6 <- cbind(id = "P-4 (descrittiva)", archetype = "A6", tissue = TIS[["A6"]], s1("K2_stardist", "A6"))
tis <- gen[gen$method != "stardist_he" & gen$phimode == "tissue", ]
roi_mean <- function(col) aggregate(tis[[col]], by = list(archetype = tis$archetype, roi = tis$roi), FUN = mean)
b4 <- roi_mean("trunc_all"); nc <- roi_mean("nc95_ev")
write.csv(b4, file.path(RES, "S1.4d_p2_roi.csv"), row.names = FALSE)
prim <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
rois <- read.csv("results/R2/R2_rois_checked.csv")
real95 <- do.call(rbind, lapply(names(prim), function(A) do.call(rbind, lapply(rois$roi_id[rois$archetype == A], function(r) {
  cc <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.4/real/%s_%s_%s.rds", A, r, prim[[A]]))$cells
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  data.frame(archetype = A, roi = r, real_nc95 = mean((cc$a_n / cc$a_t)[ev] >= 0.95)) }))))
nc <- merge(nc, real95, by = c("archetype", "roi")); names(nc)[3] <- "sim_nc95"; write.csv(nc, file.path(RES, "S1.4d_p3_roi.csv"), row.names = FALSE)
for (A in names(TIS)) { v <- b4$x[b4$archetype == A]
  P[[length(P) + 1]] <- cbind(id = "P-2", archetype = A, tissue = TIS[[A]], judge(mean(v), sd(v) / sqrt(length(v)), 0, 0.05), n_units = length(v), n_windows = NA) }
for (A in names(TIS)) { v <- nc$sim_nc95[nc$archetype == A]
  P[[length(P) + 1]] <- cbind(id = "P-3", archetype = A, tissue = TIS[[A]], judge(mean(v), sd(v) / sqrt(length(v)), 0, max(nc$real_nc95[nc$archetype == A])),
                              n_units = length(v), n_windows = NA) }
mem_gb <- { l <- readLines(file.path(OUT, "perf.log")); as.numeric(sub(".*: ", "", grep("Maximum resident", l, value = TRUE))) / 1024^2 }
pt <- judge(mean(pf$times), sd(pf$times) / sqrt(3), 0, 300, frac = FALSE)
P[[length(P) + 1]] <- cbind(id = "P-5 (tempo, s)", archetype = "S", tissue = "sezione S", pt, n_units = 3, n_windows = NA)
P[[length(P) + 1]] <- cbind(id = "P-5 (memoria, GB)", archetype = "S", tissue = "sezione S", judge(mem_gb, 0, 0, 8, frac = FALSE), n_units = 1, n_windows = NA)
pred <- do.call(rbind, P); write.csv(pred, file.path(RES, "S1.4d_predictions.csv"), row.names = FALSE)
write.csv(data.frame(times_s = pf$times, warmup_s = pf$warmup_s, mem_gb = mem_gb, env_before = paste(pf$env_before, collapse = " | "),
                     env_after = paste(pf$env_after, collapse = " | ")), file.path(RES, "S1.4d_perf.csv"), row.names = FALSE)
cp <- do.call(rbind, lapply(names(TIS), function(A) cbind(id = "CP-1", archetype = A, tissue = TIS[[A]], s1("Kc", A))))
cp$dense <- cp$archetype %in% c("A1", "A2", "A3", "A4")
write.csv(rbind(cp[, names(cp) != "dense"], d_a6), file.path(RES, "S1.4d_counterproof.csv"), row.names = FALSE)
nbad <- sum(cp$dense & cp$verdict != "smentita")
G("CP-1", "vecchio modello non smentito nei tessuti densi (deve essere < 2)", nbad, "< 2", nbad < 2)
gt <- do.call(rbind, gates); write.csv(gt, file.path(RES, "S1.4d_gates.csv"), row.names = FALSE)
print(gt, right = FALSE); print(pred[, c("id", "archetype", "est", "se", "lo2", "hi2", "int_lo", "int_hi", "verdict")], digits = 4)
print(cp[, c("id", "archetype", "est", "se", "int_lo", "int_hi", "verdict")], digits = 4)
