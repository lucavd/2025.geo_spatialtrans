# tools/S1.5d_export_presets.R - S1.5d: esporta le tabelle dei preset per tipo (.TP_TABLE, .TP_TYPES) per tools/S1.5d_size_resid_table.py.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5d_export_presets.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
source("R/04b0_tissue_presets.R")
m <- merge(.TP_TYPES[, c("preset", "cell_type", "fraction", "size_group", "rho_g", "beta_g", "sigma_g")], .TP_TABLE[, c("preset", "archetype", "density_model")], by = "preset")
write.csv(m[order(m$archetype, m$cell_type), ], "results/S1.5d/S1.5d_preset_types.csv", row.names = FALSE); cat("ok", nrow(m), "\n")
