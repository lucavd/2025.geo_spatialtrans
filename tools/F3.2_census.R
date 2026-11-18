# tools/F3.2_census.R - F3.2, cancello G-1: ogni valore che il simulatore usa per difetto nei sei tessuti, riletto dal pacchetto.
# Uso: Rscript --vanilla tools/F3.2_census.R [out.csv]   (predefinito results/F3.2/F3.2_census_values.csv)
# Righe: famiglia, funzione, tessuto, parametro, valore (scalari) oppure n e md5 dei valori (tabelle: profili, residui, quantili).
# Nessun calcolo sui dati veri: solo lettura delle costanti e dei file dei profili che il pacchetto carica.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE); out <- if (length(a)) a[1] else "results/F3.2/F3.2_census_values.csv"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
A <- paste0("A", 1:6); P <- tissue_presets(); names(P) <- A
rows <- list()
add <- function(family, fun, arch, param, value, src = NA_character_, kind = "scalar", n = 1L, md5 = NA_character_) {
  rows[[length(rows) + 1L]] <<- data.frame(family = family, fun = fun, archetype = arch, parameter = param,
                                           value = if (kind == "scalar") format(value, digits = 17) else NA_character_,
                                           kind = kind, n_values = n, values_md5 = md5, src_in_code = src, stringsAsFactors = FALSE)
}
md5v <- function(x) { f <- tempfile(); writeLines(format(unlist(x), digits = 17), f); unname(tools::md5sum(f)) }
for (ar in A) {
  sp <- tissue_preset_spec(P[[ar]])                                   # predefiniti: multi, recentred, uniform, gaussian
  s <- sp$sources
  for (i in seq_len(nrow(s))) add("cell_layer.tissue", "tissue_preset_spec()$sources", ar, s$parameter[i], s$value[i], s$source[i])
  ct <- sp$cell_types
  for (j in seq_len(nrow(ct))) for (col in c("density", "nucleus_to_eq_ratio", "nucleus_area_slope", "nucleus_area_sdlog", "nucleus_nc_max", "pale_fraction", "min_dist_um")) {
    add("cell_layer.type_catalogue", "tissue_preset_spec()$cell_types", ar, paste0(ct$cell_type[j], ":", col), ct[[col]][j],
        sp$type_sources$source[match(ct$cell_type[j], sp$type_sources$cell_type)])
  }
  co <- sp$composition
  for (j in seq_len(nrow(co))) add("cell_layer.composition", "tissue_preset_spec()$composition", ar, paste0(co$cell_type[j], ":fraction"), co$fraction[j],
                                   sp$type_sources$source[match(co$cell_type[j], sp$type_sources$cell_type)])
  cp <- sp$compartments
  for (j in seq_len(nrow(cp))) for (col in c("ell_um", "s")) add("cell_layer.compartments(option)", "tissue_preset_spec()$compartments", ar,
                                                                  paste0(cp$cell_type[j], ":", col), cp[[col]][j], cp$source[j])
  hz <- sp$heterogeneity
  for (col in setdiff(names(hz), "src")) add("cell_layer.heterogeneity", "tissue_preset_spec()$heterogeneity", ar, col, hz[[col]][1], hz$src[1])
  tp <- transcript_params(ar)
  for (nm in setdiff(names(tp), "archetype")) add("step2.transcripts", "transcript_params()", ar, nm, tp[[nm]])
  ca <- capture_params(ar)
  for (nm in setdiff(names(ca), c("archetype", "section", "parity_field", "slow"))) add("step2.capture", "capture_params()", ar, nm, ca[[nm]])
  if (!is.null(ca$parity_field)) for (nm in names(ca$parity_field)) add("step2.capture", "capture_params()$parity_field", ar, nm, ca$parity_field[[nm]])
  add("step2.depth", "visium_hd_depth_factor()", ar, "depth_factor", visium_hd_depth_factor(ar))
  pr <- load_cell_profiles(ar)
  add("step2.profiles", "load_cell_profiles()$P", ar, "P (subtype x gene)", NA, kind = "table", n = length(pr$P), md5 = md5v(round(pr$P, 15)))
  add("step2.profiles", "load_cell_profiles()$phi", ar, "phi (per gene)", NA, kind = "table", n = length(pr$phi), md5 = md5v(pr$phi))
  add("step2.profiles", "load_cell_profiles()$totals", ar, "totals (quantiles per type)", NA, kind = "table", n = prod(dim(pr$totals)), md5 = md5v(pr$totals))
  add("step2.profiles", "load_cell_profiles()$subtypes", ar, "subtype shares", NA, kind = "table", n = nrow(pr$subtypes), md5 = md5v(pr$subtypes$share))
  add("step2.profiles", "load_cell_profiles()$meta", ar, "depth_factor_S2.2 (superseded by visium_hd_depth_factor)", pr$depth_factor)
  add("step2.profiles", "load_cell_profiles()$meta", ar, "n_cells_source", pr$meta$n_cells_source)
}
rs <- .tp_resid_table()
for (ar in A) { k <- grepl(paste0("^", P[[ar]], ":"), rs[[1]]) | (if ("preset" %in% names(rs)) rs$preset == P[[ar]] else FALSE)
  add("cell_layer.size_resid(option)", "inst/extdata/size_resid.csv", ar, "empirical residuals", NA, kind = "table", n = sum(k), md5 = md5v(rs[k, ])) }
add("cell_layer.defaults", ".TP_SPATIAL_DEFAULT", "all", "spatial", NA, kind = "setting"); rows[[length(rows)]]$value <- .TP_SPATIAL_DEFAULT
add("cell_layer.defaults", ".TP_SIZE_RESID_DEFAULT", "all", "size_resid", NA, kind = "setting"); rows[[length(rows)]]$value <- .TP_SIZE_RESID_DEFAULT
add("step2.grid", "build_hd_grid() formals", "all", "margin_um", eval(formals(build_hd_grid)$margin_um))
add("step2.grid", "build_hd_grid() formals", "all", "bin_um", eval(formals(build_hd_grid)$bin_um))
add("step2.pieces", ".P27_TILE_UM", "all", "tile_um", .P27_TILE_UM)
d <- do.call(rbind, rows)
utils::write.csv(d, out, row.names = FALSE)
cat("census values:", nrow(d), "rows;", sum(d$kind == "scalar"), "scalars;", sum(d$kind == "table"), "tables ->", out, "\n")
print(table(d$family, d$archetype))
