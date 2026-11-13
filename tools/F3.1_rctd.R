# tools/F3.1_rctd.R - F3.1: RCTD (spacexr 2.2.1, commit ccc4344) in modalita' doppietto sui blocchi di 8 um di un'area, e misure.
# Riferimento: tools/F3.1_reference.py (uguale per i due lati). UMI_min = 20 per blocco e per cellula del riferimento (soglia delle
# etichette di R6); gli altri parametri ai valori predefiniti. Uso: Rscript --vanilla tools/F3.1_rctd.R <area_dir> <ref_dir> <A> <sim|real> [reuse]
# Con 'reuse' le misure si ricalcolano da rctd_blocks.csv senza rifare RCTD.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(spacexr); library(jsonlite) })
a <- commandArgs(TRUE); d <- a[1]; ref_dir <- a[2]; A <- a[3]; side <- a[4]
t0 <- proc.time()[["elapsed"]]
blk <- read.csv(file.path(d, "blocks.csv"), stringsAsFactors = FALSE)
for (v in intersect(c("nuc_type", "true_type"), names(blk))) blk[[v]][!is.na(blk[[v]]) & blk[[v]] == ""] <- NA
Xb <- as(readMM(file.path(d, "blocks_counts.mtx")), "CsparseMatrix"); g_area <- readLines(file.path(d, "genes.tsv"))
rownames(Xb) <- make.unique(g_area); colnames(Xb) <- blk$block
Xr <- as(readMM(file.path(ref_dir, sprintf("ref_%s_counts.mtx", A))), "CsparseMatrix")
rc <- read.csv(file.path(ref_dir, sprintf("ref_%s_cells.csv", A)), stringsAsFactors = FALSE)
rownames(Xr) <- make.unique(readLines(file.path(ref_dir, sprintf("ref_%s_genes.tsv", A)))); colnames(Xr) <- rc$cell_id
reuse <- length(a) >= 5 && a[5] == "reuse" && file.exists(file.path(d, "rctd_blocks.csv"))
common <- intersect(rownames(Xr), rownames(Xb))
if (reuse) {
  m <- read.csv(file.path(d, "rctd_blocks.csv"), stringsAsFactors = FALSE)
  for (v in intersect(c("nuc_type", "true_type", "spot_class", "first", "second"), names(m))) m[[v]][!is.na(m[[v]]) & m[[v]] == ""] <- NA
} else {
  ct <- factor(make.names(rc$cell_type)); names(ct) <- rc$cell_id
  type_map <- setNames(rc$cell_type, make.names(rc$cell_type))
  ref <- Reference(Xr[common, ], ct, nUMI = setNames(colSums(Xr), rc$cell_id), min_UMI = 20, n_max_cells = 10000)
  coords <- data.frame(x = blk$x, y = blk$y, row.names = blk$block)
  puck <- SpatialRNA(coords, Xb[common, ], nUMI = setNames(colSums(Xb), blk$block))
  # regola F3.1 (E-8): soglia per la stima del rumore = min(300, UMI del centesimo blocco piu' ricco); cambia solo dove meno di 100 blocchi superano 300
  ums <- min(300, sort(colSums(Xb), decreasing = TRUE)[min(100L, ncol(Xb))])
  # seme dichiarato (G-3: senza seme RCTD non e' riproducibile, perche' estrae a caso i blocchi su cui stima il rumore)
  set.seed(20261118L, kind = "Mersenne-Twister")
  rc_obj <- create.RCTD(puck, ref, max_cores = 1, UMI_min = 20, CELL_MIN_INSTANCE = 25, UMI_min_sigma = ums)
  rc_obj <- run.RCTD(rc_obj, doublet_mode = "doublet")
  rd <- rc_obj@results$results_df
  rd$block <- rownames(rd); rd$first <- type_map[as.character(rd$first_type)]; rd$second <- type_map[as.character(rd$second_type)]
  m <- merge(blk, rd[, c("block", "spot_class", "first", "second")], by = "block", all.x = TRUE)
}
run <- !is.na(m$spot_class); ok <- run & m$spot_class != "reject"
res <- list(area = basename(d), archetype = A, side = side, UMI_min_sigma = if (reuse) NA else ums, n_blocks = nrow(m), n_run = sum(run), n_genes_common = length(common),
            frac_lowumi = mean(!run), frac_reject = mean(m$spot_class[run] == "reject"),
            frac_singlet = mean(m$spot_class[ok] == "singlet"), frac_doublet_certain = mean(m$spot_class[ok] == "doublet_certain"))
hn <- ok & !is.na(m$nuc_type) & m$nuc_type %in% rc$cell_type
res$n_nuc_blocks <- sum(hn); res$agree_nucleus <- mean(m$first[hn] == m$nuc_type[hn])
# controprova: tipo del nucleo dominante del blocco spostato di 3 blocchi (24 um) in diagonale; la concordanza deve scendere
nt <- setNames(m$nuc_type, paste(m$r8, m$c8)); sh <- nt[paste(m$r8 - 3, m$c8 - 3)]
hs <- ok & !is.na(sh) & sh %in% rc$cell_type; res$agree_nucleus_shift <- mean(m$first[hs] == sh[hs]); res$n_shift_blocks <- sum(hs)
key <- paste(m$r8, m$c8); ft <- setNames(ifelse(ok, m$first, NA), key)
nb <- c(ft[paste(m$r8, m$c8 + 1)], ft[paste(m$r8 + 1, m$c8)]); self <- c(ft, ft)
pp <- !is.na(nb) & !is.na(self); res$neighbour_same <- mean(nb[pp] == self[pp]); res$n_pairs <- sum(pp)
tf <- table(factor(m$first[ok], levels = sort(unique(rc$cell_type)))); res$first_type_share <- as.list(round(tf / sum(tf), 4))
if (side == "sim") {
  tt <- ok & !is.na(m$true_type)
  res$truth_acc_first <- mean(m$first[tt] == m$true_type[tt])
  res$truth_acc_nucleus <- mean(m$nuc_type[hn & tt] == m$true_type[hn & tt])
  res$truth_acc_first_on_nuc <- mean(m$first[hn & tt] == m$true_type[hn & tt])
}
res$sec <- if (reuse && file.exists(file.path(d, "rctd.json"))) read_json(file.path(d, "rctd.json"))$sec else proc.time()[["elapsed"]] - t0
if (!reuse) write.csv(m, file.path(d, "rctd_blocks.csv"), row.names = FALSE)
write_json(res, file.path(d, "rctd.json"), auto_unbox = TRUE, digits = NA)
cat(toJSON(res[setdiff(names(res), "first_type_share")], auto_unbox = TRUE, digits = 4), "\n")
