# tools/S2.7_PH_A1.R - S2.7, anomalia A-1 (analisi post hoc): quanto la chiave unica key * (G + 1) + gene della funzione di S2.3, prima della
# correzione, fondeva geni diversi sulle aree grandi. Si prendono le terne (quadratino, cellula, gene) di un'area di 1 mm^2 del pannello
# (posizione e spostamento con la funzione corretta, strati come simulate_visium_hd()) e si rifa' il raggruppamento vecchio spostando la
# chiave alla grandezza di un'area grande: key' = key + (b - 1) * n_grande, con b estratto uniforme fra i quadratini dell'area grande.
# Uso (root del repo): Rscript --vanilla tools/S2.7_PH_A1.R <area_tag> <n_cells_grande> <n_bins_grande> <etichetta> <out.csv>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE); tag <- a[1]; nG <- as.numeric(a[2]); nbG <- as.numeric(a[3]); lab <- a[4]; OUT <- a[5]
A <- strsplit(tag, "_")[[1]][1]
L <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.6/layers/%s_U6.rds", tag)); hd <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S2.3/panel/hd/%s.rds", tag))
cp <- capture_params(A); p <- transcript_params(A); p$w2 <- cp$w2; p$ell2 <- cp$ell2
St <- simulate_stripes(hd, cp, seed = 11L); M <- St$M
prof <- simulate_cell_profiles(L$cell_df[, c("cell_id", "cell_type")], load_cell_profiles(A), seed = 12L, depth_factor = visium_hd_depth_factor(A) * M)
tx <- place_and_move_transcripts(L, hd, prof, p, seed = 13L, frame_open = TRUE, chunk_size = 2000L)$transcripts
G <- nrow(prof$counts); n <- nrow(L$cell_df)
ci <- match(tx$cell_id, L$cell_df$cell_id)
key0 <- (tx$bin_id - 1) * n + (ci - 1)
o <- order(key0, tx$gene); key0 <- key0[o]; ge <- tx$gene[o]; nn <- tx$n[o]; bi <- tx$bin_id[o]
set.seed(14L); bs <- sort(sample.int(nbG, 40L))
rows <- list()
for (b in bs) {
  key <- (b - 1 + bi - 1) * nG + (key0 %% n)                         # (quadratino, cellula) iniettivi, alla grandezza della riga b dell'area grande
  k3 <- key * (G + 1) + ge
  r <- rle(k3); e3 <- cumsum(r$lengths); rid <- rep(seq_along(e3), r$lengths)
  changed <- sum(nn[ge != ge[e3][rid]])
  rows[[length(rows) + 1L]] <- data.frame(label = lab, area = tag, b = b, k3_max = max(k3), n_triples = length(k3), n_groups = length(e3),
                                          share_triples_merged = (length(k3) - length(e3)) / length(k3), share_umi_gene_changed = changed / sum(nn))
}
R <- do.call(rbind, rows); utils::write.csv(R, OUT, row.names = FALSE)
cat(sprintf("%s %s: terne fuse in media %.4f %% (massimo %.4f %%), UMI con il gene cambiato %.4f %%\n", lab, tag, 100 * mean(R$share_triples_merged),
            100 * max(R$share_triples_merged), 100 * mean(R$share_umi_gene_changed)))
