# tools/S2.8_identity.R - S2.8, cancello G-2 e controprova CP-2: impronte delle uscite contro S2.7 e S2.5. Uso (root del repo):
#   Rscript --vanilla tools/S2.8_identity.R [results/S2.8/S2.8_G2_identity.csv]
# (a) E504 contro S2.7 p1/new (.md5: campi dei conteggi e della verita'; layers dopo aver riportato M ai nomi di S2.7 = layers_s27);
#     anche i file binari .bin .g2.bin .u8.bin .g8.bin; (b) E1008 ed EInf contro E504 (tutti i campi e i binari);
# (c) motore di S2.5 (ref) contro S2.7 p1/old (k = 2) e S2.5 panel/main_r2 (k = 0, 1): binari presenti in entrambi;
# (d) CP-2: counts_2um con ceiling_q = 1e-4 diverso da S2.7 (k = 0).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
OUTF <- commandArgs(TRUE); OUTF <- if (length(OUTF)) OUTF[1] else "results/S2.8/S2.8_G2_identity.csv"
M <- "/mnt/micron/geo_spatialtrans"; O <- file.path(M, "S2.8")
TAGS <- c("A1_r1", "A2_r1", "A3_r1", "A4_f1", "A5_r1", "A6_r1"); BINS <- c(".bin", ".g2.bin", ".u8.bin", ".g8.bin")
F27 <- c("counts_2um", "counts_8um", "bins", "bins_8um", "overlap", "overlap_8um", "origin", "origin_8um", "background", "background_8um",
         "arrived_2um", "cells", "genes")
rj <- function(f) jsonlite::fromJSON(f)
md5f <- function(f) if (file.exists(f)) unname(tools::md5sum(f)) else NA_character_
rows <- list(); add <- function(group, tag, k, field, new, ref, expect_equal = TRUE)
  rows[[length(rows) + 1L]] <<- data.frame(group = group, tag = tag, k = k, field = field, md5_new = new, md5_ref = ref,
                                            equal = !is.na(new) && !is.na(ref) && identical(new, ref), expect_equal = expect_equal)
for (t in TAGS) for (k in 0:2) {
  b <- sprintf("%s_k%d", t, k); n <- file.path(O, "e504", b); s27 <- file.path(M, "S2.7/p1/new", b)
  if (!file.exists(paste0(n, ".md5"))) { add("a_e504_vs_S2.7", t, k, "MISSING", NA, NA); next }
  hn <- rj(paste0(n, ".md5")); h7 <- rj(paste0(s27, ".md5"))
  for (fl in F27) add("a_e504_vs_S2.7", t, k, fl, hn[[fl]], h7[[fl]])
  add("a_e504_vs_S2.7", t, k, "layers(S2.7 names)", hn$layers_s27, h7$layers)
  for (e in BINS) add("a_e504_vs_S2.7", t, k, paste0("file", e), md5f(paste0(n, e)), md5f(paste0(s27, e)))
  for (cf in c("e1008", "einf")) {
    m <- file.path(O, cf, b)
    if (!file.exists(paste0(m, ".md5"))) { add(paste0("b_", cf, "_vs_e504"), t, k, "MISSING", NA, NA); next }
    hm <- rj(paste0(m, ".md5"))
    for (fl in names(hn)) add(paste0("b_", cf, "_vs_e504"), t, k, fl, hm[[fl]], hn[[fl]])
    for (e in BINS) add(paste0("b_", cf, "_vs_e504"), t, k, paste0("file", e), md5f(paste0(m, e)), md5f(paste0(n, e)))
  }
  r <- file.path(O, "ref", b); old <- if (k == 2) file.path(M, "S2.7/p1/old", b) else file.path(M, "S2.5/panel/main_r2", b)
  if (!file.exists(paste0(r, ".json"))) { add("c_ref_vs_saved", t, k, "MISSING", NA, NA); next }
  for (e in BINS) if (file.exists(paste0(old, e))) add("c_ref_vs_saved", t, k, paste0("file", e), md5f(paste0(r, e)), md5f(paste0(old, e)))
}
for (t in TAGS) {
  b <- sprintf("%s_k0", t); f <- file.path(O, "cp2", paste0(b, ".md5"))
  if (!file.exists(f)) { add("d_cp2_vs_S2.7", t, 0L, "MISSING", NA, NA, FALSE); next }
  add("d_cp2_vs_S2.7", t, 0L, "counts_2um", rj(f)$counts_2um, rj(file.path(M, "S2.7/p1/new", paste0(b, ".md5")))$counts_2um, FALSE)
}
D <- do.call(rbind, rows); D$pass <- D$equal == D$expect_equal
dir.create(dirname(OUTF), recursive = TRUE, showWarnings = FALSE); utils::write.csv(D, OUTF, row.names = FALSE)
S <- aggregate(cbind(n = 1, pass = D$pass) ~ group, data = D, FUN = sum); print(S, row.names = FALSE)
cat(sprintf("G-2: %d/%d confronti come attesi (gruppi a-c); CP-2: %d/6 diversi\n", sum(D$pass[D$group != "d_cp2_vs_S2.7"]), sum(D$group != "d_cp2_vs_S2.7"),
            sum(D$pass[D$group == "d_cp2_vs_S2.7"])))
