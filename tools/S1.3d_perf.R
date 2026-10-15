# tools/S1.3d_perf.R — S1.3d: tempo e memoria su S (pattern C-perf, 449 536 celle) per candidato x operazione, processo separato.
# Uso: /usr/bin/time -v Rscript --vanilla tools/S1.3d_perf.R <cand> <op = O1 | O2s2 | O2s8> <rep> [outdir]
# Stampa PERF,<cand>,<op>,<rep>,<t_prepare>,<t_op>; con rep = 1 salva il confronto con la verità in res/S_cperf_<cand>_<op>.rds.
# F4: t_prepare include f4_build() (costruzione della tassellazione dai generatori).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
a <- commandArgs(TRUE); CAND <- a[1]; OP <- a[2]; REP <- as.integer(a[3]); OUTD <- if (length(a) > 3) a[4] else "/mnt/micron/geo_spatialtrans/S1.3d"
source("tools/S1.3d_common.R"); source("tools/S1.3d_candidates.R"); source("tools/S1.3d_f4.R"); source("tools/S1.3d_compare.R")
suppressPackageStartupMessages(library(jsonlite))
GEOM <- file.path(OUTD, "geom"); TRUTH <- file.path(OUTD, "truth"); RES <- file.path(OUTD, "res"); key <- "S_cperf"
meta <- fromJSON(file.path(GEOM, paste0(key, "_meta.json")), simplifyVector = TRUE)
terr <- df_to_sfc(as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_cells.parquet")))), meta$n)
reg <- df_to_sfc(transform(as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_region.parquet")))), cell = 1L), 1L)
pairs <- as.matrix(as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_pairs.parquet")))))
pts <- as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_pts.parquet"))))
t0 <- proc.time()[[3]]; f4 <- NULL; f4info <- NULL
if (CAND == "F4") { b <- f4_build(pts$x, pts$y, reg); f4 <- b$terr; f4info <- b$info[setdiff(names(b$info), c("lost", "gained"))]
  f4info$same_lost <- identical(as.integer(b$info$lost), as.integer(meta$n_lost)); f4info$same_gained <- identical(as.integer(b$info$gained), as.integer(meta$n_gained))
  a0 <- as.numeric(st_area(terr)); a4 <- as.numeric(st_area(f4)); f4info$c5_n_over <- sum(abs(a4 - a0) > G_GRID * meta$P); f4info$c5_max_ratio <- max(abs(a4 - a0) / (G_GRID * meta$P)) }
C <- cand_prepare(CAND, terr, reg, f4 = f4); t_prep <- proc.time()[[3]] - t0
if (CAND == "F2") { f2info <- list(n_invalid_rounded = length(C$info$invalid_idx), still_invalid = C$info$still_invalid)
  if (identical(Sys.getenv("S13D_F2INFO_ONLY"), "1")) { saveRDS(f2info, file.path(RES, "S_cperf_F2_validity.rds")); print(f2info); quit(status = 0) } }
t1 <- proc.time()[[3]]
if (OP == "O1") o <- o1_run(C, pairs) else { s <- if (OP == "O2s2") 2 else 8; o <- o2_run(C, meta$gs, s) }
t_op <- proc.time()[[3]] - t1
cat(sprintf("PERF,%s,%s,%d,%.2f,%.2f\n", CAND, OP, REP, t_prep, t_op))
if (REP == 1) {
  if (OP == "O1") {
    tr <- fromJSON(file.path(TRUTH, paste0(key, "_o1.json")))$o1
    tp <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o1_pairs.parquet")))); tc <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o1_cells.parquet"))))
    stopifnot(identical(as.integer(tp$i), as.integer(pairs[, 1])))
    cm <- cmp_o1(o, tr, tp, tc, meta$P, meta$Preg)
  } else {
    tt <- as.data.frame(read_parquet(file.path(TRUTH, sprintf("%s_o2_s%d.parquet", key, s))))
    cm <- cmp_o2(o$pieces, tt, s, CAND)
    if (s == 2) { t8 <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o2_s8.parquet"))))
      k8 <- paste(tt$i, tt$ix %/% 4, tt$iy %/% 4); tau3 <- tapply(G_GRID * tt$perim, k8, sum)[paste(t8$i, t8$ix, t8$iy)]
      cm$o3 <- cmp_o2(o3_aggregate(o$pieces), t8, 8, CAND, tau_truth = as.numeric(tau3)) }
  }
  if (CAND == "F4" && OP == "O1") write_parquet(geom_to_df(f4), file.path(GEOM, paste0(key, "_f4cells.parquet")))
  saveRDS(list(cand = CAND, op = OP, t_prepare = t_prep, t_op = t_op, cmp = cm, f4info = f4info), file.path(RES, sprintf("%s_%s_%s.rds", key, CAND, OP)))
}
