# tools/F3.1_export_sim.R - F3.1: dall'uscita di simulate_visium_hd_pieces() al formato comune dei due lati (un'area).
# Uso: Rscript --vanilla tools/F3.1_export_sim.R <rds> <out_dir> <side_um>
# Scrive in out_dir: X.parquet (terne quadratino, gene, UMI; indici 0-based nell'ordine di bins/genes), bins.parquet
# (bin_id, row, col, x_um, y_um, in_tissue, label = cellula il cui nucleo copre >= meta' del quadratino, 0 altrimenti),
# cells.parquet (label, cell_type), genes.tsv, profiles.parquet (profili veri generati: label, gene, n), origin.parquet (bin_id, label, n: trascritti catturati per coppia; senza fondo).
# Quadratini: centro dentro il quadrato [0, side] x [0, side] dell'area (stesso riquadro del lato vero).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(arrow) })
a <- commandArgs(TRUE); o <- readRDS(a[1]); out <- a[2]; side <- as.numeric(a[3])
dir.create(out, recursive = TRUE, showWarnings = FALSE)
b <- o$bins
keep <- b$x >= 0 & b$x <= side & b$y >= 0 & b$y <= side
b <- b[keep, ]
ov <- o$truth$overlap
ov <- ov[ov$frac_nucleus >= 0.5 & ov$bin_id %in% b$bin_id, c("bin_id", "cell_id", "frac_nucleus")]
ov <- ov[order(ov$bin_id, -ov$frac_nucleus), ]; ov <- ov[!duplicated(ov$bin_id), ]
lab <- integer(nrow(b)); lab[match(ov$bin_id, b$bin_id)] <- as.integer(ov$cell_id)
bins <- data.frame(bin_id = as.integer(b$bin_id), row = as.integer(b$i), col = as.integer(b$j), x_um = b$x, y_um = b$y,
                   in_tissue = b$in_tissue, label = lab)
write_parquet(bins, file.path(out, "bins.parquet"))
X <- o$counts_2um[, b$bin_id, drop = FALSE]          # geni x quadratini (colonne = bin_id, verificato sotto)
stopifnot(ncol(o$counts_2um) == nrow(o$bins), all(o$bins$bin_id == seq_len(nrow(o$bins))))
T <- as(X, "TsparseMatrix")
write_parquet(data.frame(bin = T@j, gene = T@i, n = as.integer(T@x)), file.path(out, "X.parquet"))
writeLines(rownames(o$counts_2um), file.path(out, "genes.tsv"))
cl <- o$truth$cells
# profili veri generati (geni x cellule, prima di spostamento e cattura) delle cellule con almeno un quadratino di nucleo nel riquadro
cc <- o$truth$cell_counts; ids <- sort(unique(lab[lab > 0])); Pm <- as(cc[, match(ids, cl$cell_id), drop = FALSE], "TsparseMatrix")
write_parquet(data.frame(label = ids[Pm@j + 1L], gene = Pm@i, n = as.integer(round(Pm@x))), file.path(out, "profiles.parquet"))
write_parquet(data.frame(label = as.integer(cl$cell_id), cell_type = cl$cell_type, subtype = cl$subtype), file.path(out, "cells.parquet"))
og <- o$truth$origin; og <- og[og$bin_id %in% b$bin_id, ]
write_parquet(data.frame(bin_id = as.integer(og$bin_id), label = as.integer(og$cell_id), n = as.integer(og$n)), file.path(out, "origin.parquet"))
bg <- o$truth$background; bg <- bg[bg$bin_id %in% b$bin_id, ]
write_parquet(as.data.frame(bg), file.path(out, "background.parquet"))
chk <- c(n_bins = nrow(bins), umi = sum(T@x), umi_origin = sum(og$n), umi_bg = sum(bg[[ncol(bg)]]), n_lab_bins = sum(lab > 0),
         n_cells = nrow(cl))
jsonlite::write_json(c(as.list(chk), list(tag = o$f31$tag, k = o$f31$k, seed = o$f31$seed, archetype = o$f31$archetype, side_um = side,
                     gene_head = head(rownames(o$counts_2um), 3), bg_cols = names(bg))), file.path(out, "export.json"), auto_unbox = TRUE)
cat(o$f31$tag, "bins", nrow(bins), "umi", sum(T@x), "origin+bg", sum(og$n) + sum(bg[[ncol(bg)]]), "lab bins", sum(lab > 0), "\n")
