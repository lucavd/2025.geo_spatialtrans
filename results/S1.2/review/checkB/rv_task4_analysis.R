# rv_task4_analysis.R — revisione avversariale S1.2 / checkB, compito 4: analisi della sensibilita' alla banda.
# Legge sens/*.rds (rv_task4_sens.R). Codice base R scritto da zero (nessuna funzione del progetto).
# (a) riproducibilita': variante x1 vs /mnt/micron/.../sim/*.rds e real_pcf.rds del progetto;
# (b) CP-1, CP-2 (e CP-3) ricalcolati per ogni variante di banda: x1, x0.5, x2, own (Stoyan per pattern).
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task4_analysis.R
SO <- "results/S1.2/review/checkB/sens"; RV <- "results/S1.2/review/checkB"
R_GRID <- seq(0, 30, by = 0.5); I <- R_GRID > 0
fs <- list.files(SO, pattern = "\\.rds$", full.names = TRUE)
S <- do.call(rbind, lapply(fs, readRDS))
cat("unita' lette:", length(fs), " righe:", nrow(S), "\n")
trapD <- function(a, b, idx = I) { f <- (a - b)^2; x <- R_GRID[idx]; f <- f[idx]; sum((head(f, -1) + tail(f, -1)) / 2 * diff(x)) }
key <- function(...) paste(..., sep = "|")
S$k <- key(S$archetype, S$roi_id, S$source, S$variant)
G <- split(S$g, S$k)            # ordine di r preservato per costruzione
meta <- unique(S[, c("archetype", "roi_id", "source", "variant", "d", "feasible", "reps", "chk_maxabs", "bw0")])
# (a) riproducibilita'
OUT <- "/mnt/micron/geo_spatialtrans/S1.2/checkB"
ps <- do.call(rbind, lapply(list.files(file.path(OUT, "sim"), full.names = TRUE), readRDS))
gc_ <- grep("^g_", names(ps), value = TRUE)
rep_rows <- lapply(seq_len(nrow(ps)), function(i) { k <- key(ps$archetype[i], ps$roi_id[i], ps$model[i], "x1")
  m <- meta[meta$archetype == ps$archetype[i] & meta$roi_id == ps$roi_id[i] & meta$source == ps$model[i] & meta$variant == "x1", ]
  if (!k %in% names(G)) return(data.frame(archetype = ps$archetype[i], roi_id = ps$roi_id[i], model = ps$model[i], present = FALSE, maxabs = NA, same_feasible = NA, same_reps = NA))
  data.frame(archetype = ps$archetype[i], roi_id = ps$roi_id[i], model = ps$model[i], present = TRUE,
             maxabs = if (ps$reps[i] > 0) max(abs(unlist(ps[i, gc_]) - G[[k]])[I]) else NA, same_feasible = ps$feasible[i] == m$feasible, same_reps = ps$reps[i] == m$reps) })
rep_df <- do.call(rbind, rep_rows)
real_ref <- readRDS(file.path(OUT, "real_pcf.rds"))
rr <- do.call(rbind, lapply(real_ref, function(z) data.frame(archetype = z$archetype, roi_id = z$roi_id,
  maxabs_primary = max(abs(z$primary$g - G[[key(z$archetype, z$roi_id, "real_primary", "x1")]])[I]),
  maxabs_secondary = max(abs(z$secondary$g - G[[key(z$archetype, z$roi_id, "real_secondary", "x1")]])[I]))))
write.csv(rep_df, file.path(RV, "rv_task4_repro_sim.csv"), row.names = FALSE); write.csv(rr, file.path(RV, "rv_task4_repro_real.csv"), row.names = FALSE)
cat(sprintf("riproducibilita' sim: %d/%d presenti, max|dg| %.2e, feasible uguali %d/%d, reps uguali %d/%d; reale max|dg| %.2e / %.2e; chk vs spatstat max %.2e\n",
    sum(rep_df$present), nrow(rep_df), max(rep_df$maxabs, na.rm = TRUE), sum(rep_df$same_feasible, na.rm = TRUE), sum(rep_df$present),
    sum(rep_df$same_reps, na.rm = TRUE), sum(rep_df$present), max(rr$maxabs_primary), max(rr$maxabs_secondary), max(meta$chk_maxabs, na.rm = TRUE)))
# (b) controprove per variante
ARCH <- data.frame(archetype = paste0("A", 1:6), tot = c(12419, 8924, 3100, 28096, 1185, 988), d_nuc = c(4.1, 6.3, 6.3, 4.2, 11.3, 4.8))
ARCH$d_rule <- 2 / 3 * sqrt(1e6 / (pi * ARCH$tot))
rows <- list(); cp2c <- list()
for (v in c("x1", "x0.5", "x2", "own", "x1_from_r0")) {
  vv <- if (v == "x1_from_r0") "x1" else v; idx <- if (v == "x1_from_r0") rep(TRUE, length(R_GRID)) else I
  Dm <- function(A, roi, src) { k <- key(A, roi, src, vv); g <- G[[k]]; if (is.null(g) || all(is.na(g))) NA else trapD(g, G[[key(A, roi, "real_primary", vv)]], idx) }
  just <- 0; p3 <- 0
  for (A in ARCH$archetype) {
    a <- ARCH[ARCH$archetype == A, ]; rois <- sort(unique(meta$roi_id[meta$archetype == A]))
    Dpd <- sapply(rois, Dm, A = A, src = "PD"); Dcs <- sapply(rois, Dm, A = A, src = "CSR"); Dnu <- sapply(rois, Dm, A = A, src = "NUC")
    red <- 1 - Dpd / Dcs; j <- sum(Dpd < Dcs) >= 4 && median(red) >= 0.25; just <- just + j
    dg <- seq(0, 15, 0.5)
    Dmed <- sapply(dg, function(d) { src <- sprintf("G%04.1f", d)
      ok <- all(sapply(rois, function(roi) { m <- meta[meta$archetype == A & meta$roi_id == roi & meta$source == src & meta$variant == vv, ]; nrow(m) == 1 && m$feasible && m$reps == 20 }))
      if (ok) median(sapply(rois, Dm, A = A, src = src)) else NA })
    cp2c[[length(cp2c) + 1]] <- data.frame(variant = v, archetype = A, d = dg, D_med = Dmed)
    ds <- dg[which.min(Dmed)]; Ds <- min(Dmed, na.rm = TRUE); ratio <- median(Dpd) / Ds
    Dseg <- sapply(rois, function(roi) trapD(G[[key(A, roi, "real_secondary", vv)]], G[[key(A, roi, "real_primary", vv)]], idx))
    c3 <- median(Dseg) / median(Dcs); p3 <- p3 + (c3 < 0.5)
    rows[[length(rows) + 1]] <- data.frame(variant = v, archetype = A, cp1_wins = sum(Dpd < Dcs), cp1_red_med = median(red), cp1 = ifelse(j, "PASS", "FAIL"),
      d_star = ds, D_star = Ds, d_rule = a$d_rule, ratio_rule = ratio, cp2 = ifelse(ratio <= 1.25, "PASS", "FAIL"), ratio_nuc = median(Dnu) / Ds,
      d_max_feasible = max(dg[!is.na(Dmed)]), cp3_ratio = c3, cp3 = ifelse(c3 < 0.5, "PASS", "FAIL"))
  }
  rows[[length(rows) + 1]] <- data.frame(variant = v, archetype = "tutti", cp1_wins = NA, cp1_red_med = NA, cp1 = sprintf("%d/6 %s", just, ifelse(just >= 4, "PASS", "FAIL")),
      d_star = NA, D_star = NA, d_rule = NA, ratio_rule = NA, cp2 = NA, ratio_nuc = NA, d_max_feasible = NA, cp3_ratio = NA, cp3 = sprintf("%d/6", p3))
}
res <- do.call(rbind, rows); write.csv(res, file.path(RV, "rv_task4_sens_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, cp2c), file.path(RV, "rv_task4_cp2_curves.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
