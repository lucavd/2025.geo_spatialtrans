# tools/sonda/S5_diag_rctd.R - S5 (docs/sonda/REGOLA_S5.md), diagnosi dopo il pilota, nessuna scelta di analisi: per ogni esecuzione della lista
# ricostruisce gli oggetti di tools/F3.1_rctd.R (invariato: stesso riferimento L, UMI_min = 20, CELL_MIN_INSTANCE = 25, UMI_min_sigma come E-8,
# seme 20261118) e chiama solo create.RCTD (non run.RCTD), che contiene il controllo di spacexr get_de_genes: servono almeno il 10 % dei geni
# del riferimento con >= 3 conteggi sui blocchi con >= UMI_min UMI. Riporta blocchi totali, blocchi sopra soglia, geni con >= 3 conteggi, soglia
# e l'esito (ok o messaggio d'errore); blocchi sopra soglia = UMI sui geni comuni fra UMI_min e UMI_max (come restrict_counts). Uso: Rscript --vanilla tools/sonda/S5_diag_rctd.R <lista.tsv> <out.csv>   (lista: arm unit run_dir)
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(spacexr) })
a <- commandArgs(TRUE); L <- read.delim(a[1], header = FALSE, stringsAsFactors = FALSE); ref_dir <- "/mnt/micron/geo_spatialtrans/S3_sonda/ref"; A <- "L"
Xr <- as(readMM(file.path(ref_dir, sprintf("ref_%s_counts.mtx", A))), "CsparseMatrix"); rc <- read.csv(file.path(ref_dir, sprintf("ref_%s_cells.csv", A)), stringsAsFactors = FALSE)
rownames(Xr) <- make.unique(readLines(file.path(ref_dir, sprintf("ref_%s_genes.tsv", A)))); colnames(Xr) <- rc$cell_id
fm <- formals(create.RCTD); out <- list()
for (i in seq_len(nrow(L))) {
  d <- L$V3[i]
  if (!file.exists(file.path(d, "blocks.csv"))) { out[[i]] <- data.frame(arm = L$V1[i], unit = L$V2[i], n_blocks = NA, umi_total = NA, n_blocks_umi20 = NA, max_block_umi = NA,
      n_genes_common = NA, n_genes_ge3 = NA, threshold_genes = NA, create_ok = FALSE, message = "no blocks.csv (bin2cell failed)"); next }
  blk <- read.csv(file.path(d, "blocks.csv"), stringsAsFactors = FALSE)
  Xb <- as(readMM(file.path(d, "blocks_counts.mtx")), "CsparseMatrix"); rownames(Xb) <- make.unique(readLines(file.path(d, "genes.tsv"))); colnames(Xb) <- blk$block
  common <- intersect(rownames(Xr), rownames(Xb)); ct <- factor(make.names(rc$cell_type)); names(ct) <- rc$cell_id
  ref <- suppressWarnings(Reference(Xr[common, ], ct, nUMI = setNames(colSums(Xr), rc$cell_id), min_UMI = 20, n_max_cells = 10000))
  puck <- SpatialRNA(data.frame(x = blk$x, y = blk$y, row.names = blk$block), Xb[common, ], nUMI = setNames(colSums(Xb), blk$block))
  ums <- min(300, sort(colSums(Xb), decreasing = TRUE)[min(100L, ncol(Xb))])
  keepb <- colSums(Xb[common, ]) >= 20 & colSums(Xb[common, ]) <= fm$UMI_max
  nb20 <- sum(keepb); ng3 <- if (nb20 > 0) sum(rowSums(Xb[common, keepb, drop = FALSE]) >= 3) else 0
  set.seed(20261118L, kind = "Mersenne-Twister")
  msg <- tryCatch({ suppressMessages(suppressWarnings(create.RCTD(puck, ref, max_cores = 1, UMI_min = 20, CELL_MIN_INSTANCE = 25, UMI_min_sigma = ums))); "ok" },
                  error = function(e) conditionMessage(e))
  out[[i]] <- data.frame(arm = L$V1[i], unit = L$V2[i], n_blocks = ncol(Xb), umi_total = sum(Xb), n_blocks_umi20 = nb20, max_block_umi = max(colSums(Xb)),
                         n_genes_common = length(common), n_genes_ge3 = ng3, threshold_genes = 0.1 * length(common), create_ok = msg == "ok", message = substr(msg, 1, 120))
  cat(L$V1[i], L$V2[i], nb20, ng3, msg == "ok", "\n")
}
write.csv(do.call(rbind, out), a[2], row.names = FALSE)
