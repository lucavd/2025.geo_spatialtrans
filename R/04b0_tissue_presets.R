# R/04b0_tissue_presets.R - Step 1, S1.5: i sei preset dei tessuti (pannello A1-A6), un tipo per tessuto.
# Ogni valore ha la sua riga di docs/BIO_REFERENCES.md (colonna src_*); il test di S1.5 li rilegge dal file.
# Design §4.6 (v1.1); file dello step docs/steps/S1.5.md (definizioni operative, scelte decise dai dati).
# Densita' = totale del conteggio manuale (R2b, decisione di Luca del 2026-10-07); distanza minima = d* di S1.2 (DEC-009);
# modello di taglia M2 (S1.4, R3b/DEC-006 per linfonodo e cuore); tetto phi (S1.4c, R3b); intervallo C8 (S1.4, R3b).

.TP_TABLE <- data.frame(
  preset          = c("epithelial", "tumor_microenv", "stromal_rich", "lymphoid", "brain_cortex", "columnar_muscle"),
  archetype       = paste0("A", 1:6),
  tissue          = c("intestino di topo", "tumore del colon", "stroma del tumore del colon", "linfonodo",
                      "cervello di topo", "cuore"),
  density         = c(12419, 8924, 3100, 28096, 1185, 988),       # nuclei/mm², totale manuale
  src_density     = c("B-043", "B-044", "B-045", "B-046", "B-047", "B-048"),
  density_evident = c(10291, 7376, 2311, 26844, 953, 906),         # nuclei/mm², evidenti (stesse righe)
  min_dist_um     = c(2.5, 3.5, 3.5, 4, 7, 2.5),
  src_min_dist    = c("B-049", "B-050", "B-051", "B-052", "B-053", "B-054"),
  rho             = c(0.380, 0.510, 0.374, 0.679, 0.282, 0.154),
  src_rho         = c("B-091", "B-094", "B-097", "B-139", "B-103", "B-150"),
  beta            = c(0.114, 0.320, 0.147, 0.455, 0.621, 0.122),
  src_beta        = c("B-092", "B-095", "B-098", "B-140", "B-104", "B-151"),
  sigma           = c(0.475, 0.547, 0.553, 0.210, 0.813, 0.490),
  src_sigma       = c("B-093", "B-096", "B-099", "B-141", "B-105", "B-152"),
  # densita' con cui e' stato stimato il modello di taglia (1e6 / A_ref): results/S1.4/S1.4_calibration.csv
  # (A1, A2, A3, A5) e results/R3b/R3b_refs.csv, riga density_mm2 della maschera unica (A4, A6)
  density_model   = c(10591.0, 7723.3, 4278.9, 27897.5, 1039.7, 1062.4),
  src_density_model = c("S1.4_calibration.csv", "S1.4_calibration.csv", "S1.4_calibration.csv", "R3b_refs.csv",
                        "S1.4_calibration.csv", "R3b_refs.csv"),
  phi             = c(0.757, 1, 0.947, 0.791, 0.474, 0.149),
  src_phi         = c("B-121", "B-122", "B-123", "B-143", "B-125", "B-154"),
  c8_lo           = c(2.81, 3.30, 3.60, 2.59, 6.95, 8.84),
  c8_hi           = c(7.71, 9.72, 14.01, 4.39, 26.88, 24.01),
  src_c8          = c("B-109", "B-110", "B-111", "B-142", "B-113", "B-153"),
  stringsAsFactors = FALSE
)

#' Nomi dei preset disponibili
#' @export
tissue_presets <- function() .TP_TABLE$preset

#' Specifica di un preset di tessuto
#'
#' Un solo tipo cellulare per tessuto (`<preset>_all`, perimetro di S1.5, DEC-008): tutti i cluster ricevono quel tipo.
#' Conversione del modello di taglia alla densita' del preset (scelta 2 di S1.5): `"preserve_area"` usa
#' rho' = rho * (density / density_model)^((1 - beta) / 2), che conserva la relazione misurata fra area del nucleo e area
#' del territorio; `"none"` usa rho com'e'.
#'
#' @param preset nome del preset (vedi `tissue_presets()`).
#' @param size_conversion `"preserve_area"` (predefinito) o `"none"`.
#' @return lista con `cell_types` (catalogo di una riga), `c8_interval` (µm), `archetype`, `tissue`, `sources`
#'   (data.frame parametro, valore, fonte), `size_conversion`.
#' @export
tissue_preset_spec <- function(preset, size_conversion = c("preserve_area", "none")) {
  size_conversion <- match.arg(size_conversion)
  if (!is.character(preset) || length(preset) != 1L || !preset %in% .TP_TABLE$preset)
    stop("tissue_preset_spec(): preset sconosciuto; disponibili: ", paste(.TP_TABLE$preset, collapse = ", "), call. = FALSE)
  tb <- .TP_TABLE
  if (Sys.getenv("S15_MUTANT", "none") == "M6") tb$phi[c(4, 6)] <- tb$phi[c(6, 4)]   # mutante: phi di linfonodo e cuore scambiati
  p <- tb[tb$preset == preset, ]
  rho <- if (size_conversion == "preserve_area") p$rho * (p$density / p$density_model)^((1 - p$beta) / 2) else p$rho
  ct <- data.frame(cell_type = paste0(preset, "_all"), density = p$density, min_dist_um = p$min_dist_um,
                   nucleus_shape = "circle", nucleus_radius = NA_real_, nucleus_to_eq_ratio = rho,
                   nucleus_area_slope = p$beta, nucleus_area_sdlog = p$sigma, nucleus_aspect = 1,
                   nucleus_orient = "random", nucleus_nc_max = p$phi,
                   pale_fraction = 1 - p$density_evident / p$density, stringsAsFactors = FALSE)
  src <- data.frame(parameter = c("density", "density_evident", "min_dist_um", "rho", "beta", "sigma", "density_model",
                                  "phi", "c8_lo", "c8_hi"),
                    value = c(p$density, p$density_evident, p$min_dist_um, p$rho, p$beta, p$sigma, p$density_model,
                              p$phi, p$c8_lo, p$c8_hi),
                    source = c(p$src_density, p$src_density, p$src_min_dist, p$src_rho, p$src_beta, p$src_sigma,
                               p$src_density_model, p$src_phi, p$src_c8, p$src_c8), stringsAsFactors = FALSE)
  list(cell_types = ct, c8_interval = c(p$c8_lo, p$c8_hi), archetype = p$archetype, tissue = p$tissue,
       sources = src, size_conversion = size_conversion, preset = preset)
}
