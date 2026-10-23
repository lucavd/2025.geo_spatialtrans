# tools/S1.5_single_preset.R - dopo S1.5b i preset per tipo sono il comportamento predefinito di tissue_preset_spec() e place_cells().
# Gli script e i test di S1.5 lo caricano dopo R/04b*.R per riprodurre S1.5 con il preset a un tipo (types = "single"), senza cambiare
# nient'altro: con types = "single" il codice e' quello di S1.5 (G-3 di S1.5b: md5 identico in 300/300 unita').
.tps_s15b <- tissue_preset_spec
tissue_preset_spec <- function(preset, size_conversion = c("preserve_area", "none"), types = "single", ...) .tps_s15b(preset, size_conversion, types = types, ...)
.pc_s15b <- place_cells
place_cells <- function(..., types = "single") .pc_s15b(..., types = types)
