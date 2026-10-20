# tools/S1.4d_g4a_cells.R — S1.4d, dopo la revisione avversariale (rilievo F-04): G-4a «stesse celle al tetto» cella per cella
# su tutte le 300 unita' x seme (phi del fold), contro i nuclei K-2 del prototipo S1.4c; e md5 del seme 1 = banco (il codice
# corretto dopo la revisione da' le stesse uscite). Uso (root del repo): Rscript --vanilla tools/S1.4d_g4a_cells.R
L <- readLines("tools/S1.4d_bench.R"); k <- grep('^if \\(STAGE == "gen"\\)', L); args <- c("none", "/mnt/micron/geo_spatialtrans/S1.4d")
eval(parse(text = L[1:(k - 1)]))
one <- function(A, k) {
  method <- ARCH$primary[ARCH$archetype == A]
  o14 <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, k))); pc <- readRDS(file.path(S14C, "gen", sprintf("%s_%s_fold%d.rds", A, method, k)))
  gb <- readRDS(file.path(OUT, "gen", sprintf("%s_%s_fold%d_fold.rds", A, method, k)))
  te <- load_real(A, o14$roi_test, method); inp <- inputs_of(te); ct <- ct_m2(o14$fits$M2, pc$phi)
  do.call(rbind, lapply(seq_len(NSEED), function(s) {
    o <- E$derive_cell_geometry(inp$cen, inp$tess, ct, random_seed = o14$seed + 100000 * s, verbose = FALSE); cd <- o$cell_df
    gp <- geos_read_wkb(unclass(read_parquet(file.path(S14C, "gen_in", sprintf("%s_%s_fold%d_s%02d_K2.parquet", A, o14$roi_test, k, s)))$wkb))
    tr_p <- geos_area(gp) < (1 - 1e-6) * cd$nucleus_area_target; tr_n <- cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target
    data.frame(archetype = A, fold = k, seed_idx = s, n = nrow(cd), n_trunc_new = sum(tr_n), n_trunc_proto = sum(tr_p), n_capped_new = sum(cd$nucleus_capped),
               diff_trunc = sum(tr_n != tr_p), diff_capped_vs_proto = sum(cd$nucleus_capped != tr_p),
               md5_same_as_bench = if (s == 1) identical(wkb_md5(o), gb$md5[gb$seed_idx == 1]) else NA)
  }))
}
jobs <- expand.grid(A = ARCH$archetype, k = 1:5, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(i) one(jobs$A[i], jobs$k[i]), mc.cores = 4, mc.preschedule = FALSE)
bad <- vapply(res, inherits, TRUE, "try-error"); if (any(bad)) { print(res[bad]); quit(status = 1) }
x <- do.call(rbind, res); write.csv(x, "results/S1.4d/S1.4d_g4a_cells.csv", row.names = FALSE)
cat(sprintf("G4A_CELLS unita' x seme %d, celle %d, celle troncate diverse %d, al tetto contro troncate del prototipo diverse %d, md5 seme 1 = banco %d/%d\n",
            nrow(x), sum(x$n), sum(x$diff_trunc), sum(x$diff_capped_vs_proto), sum(x$md5_same_as_bench, na.rm = TRUE), sum(!is.na(x$md5_same_as_bench))))
