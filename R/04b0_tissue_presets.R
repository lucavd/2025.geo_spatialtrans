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

# ---- S1.5b: tipi cellulari nei preset (docs/steps/S1.5b.md; DEC-016, DEC-017, DEC-018; riferimenti B-186...B-208) ----
# Per tipo: frazione nel preset (R6, criterio corretto, fra i nuclei con un tipo; cervello corretto per i nuclei persi da Cellpose),
# gruppo di taglia e modello del gruppo (rho_g, beta_g, sigma_g con area di riferimento 1e6 / density_model del tessuto),
# territorio mediano reale del tipo (m_t) e territorio mediano simulato del tessuto in S1.5 (m_sim). Generata da
# tools/S1.5b_calibration.py (results/S1.5b/S1.5b_type_calibration.csv), 6 cifre significative.
.TP_TYPES <- data.frame(
  preset              = c("epithelial", "epithelial", "epithelial", "tumor_microenv", "tumor_microenv", "stromal_rich", "stromal_rich", "stromal_rich", "stromal_rich", "lymphoid", "lymphoid", "lymphoid", "lymphoid", "lymphoid", "lymphoid", "brain_cortex", "brain_cortex", "brain_cortex", "brain_cortex", "brain_cortex", "columnar_muscle", "columnar_muscle", "columnar_muscle"),
  cell_type           = c("epiteliale", "immunitario", "stromale e muscolo liscio", "epiteliale", "fibroblasti", "endoteliale", "epiteliale", "fibroblasti", "immunitario", "endoteliali", "linfociti B", "linfociti T", "mieloidi", "plasmacellule", "stromali", "astrociti", "microglia", "oligodendrociti", "vascolari", "neuroni", "cardiomiociti", "fibroblasti", "immunitari"),
  fraction            = c(0.685091, 0.10098, 0.213929, 0.776428, 0.223572, 0.0813992, 0.40709, 0.325079, 0.186432, 0.0152395, 0.494632, 0.416747, 0.0202196, 0.0319078, 0.0212537, 0.0610637, 0.0313428, 0.168381, 0.0396766, 0.699536, 0.757619, 0.12759, 0.114791),
  fraction_follicle   = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.0144299, 0.608699, 0.298286, 0.0240499, 0.0374975, 0.0170381, NA, NA, NA, NA, NA, NA, NA, NA),
  fraction_paracortex = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.0157831, 0.418048, 0.496282, 0.017648, 0.0281549, 0.0240841, NA, NA, NA, NA, NA, NA, NA, NA),
  size_group          = c(1, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1),
  rho_g               = c(0.398694, 0.341871, 0.341871, 0.545033, 0.43163, 0.383156, 0.383156, 0.383156, 0.338571, 0.679, 0.679, 0.679, 0.679, 0.679, 0.679, 0.19996, 0.19996, 0.19996, 0.19996, 0.324497, 0.154, 0.154, 0.154),
  beta_g              = c(0.159103, -0.00344299, -0.00344299, 0.406403, 0.322068, 0.157005, 0.157005, 0.157005, 0.160364, 0.455, 0.455, 0.455, 0.455, 0.455, 0.455, 0.588567, 0.588567, 0.588567, 0.588567, 0.20495, 0.122, 0.122, 0.122),
  sigma_g             = c(0.458977, 0.460479, 0.460479, 0.504294, 0.510044, 0.549589, 0.549589, 0.549589, 0.518084, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.867332, 0.867332, 0.867332, 0.867332, 0.478454, 0.49, 0.49, 0.49),
  m_t                 = c(63.4727, 45.7939, 67.868, 83.0454, 132.922, 164.501, 87.4305, 301.76, 186.516, 32.2683, 37.4191, 32.7613, 39.694, 34.4772, 33.0946, 740.313, 1099.8, 370.648, 1016.3, 899.307, 790.48, 618.319, 469.06),
  m_sim               = c(73.0345, 73.0345, 73.0345, 103.169, 103.169, 293.228, 293.228, 293.228, 293.228, 34.8323, 34.8323, 34.8323, 34.8323, 34.8323, 34.8323, 746.659, 746.659, 746.659, 746.659, 746.659, 870.538, 870.538, 870.538),
  src                 = c("B-186", "B-187", "B-188", "B-189", "B-190", "B-191", "B-192", "B-193", "B-194", "B-195", "B-196", "B-197", "B-198", "B-199", "B-200", "B-201", "B-202", "B-203", "B-204", "B-205", "B-206", "B-207", "B-208"),
  stringsAsFactors = FALSE
)

# ---- S1.5c: compartimenti per tipo (docs/steps/S1.5c.md; riferimenti B-209...B-231) ----
# Per tipo: lunghezza di correlazione ell_um e ampiezza s del campo latente (calibrazione su tutte le 5 ROI, tools/S1.5c_calibration.py;
# tipi minori s = 0, senza campo) e area nucleare mediana reale del tipo (tools/S1.5c_type_nuc_area.py), per la distanza minima per nucleo.
.TP_COMPART <- data.frame(
  preset = c("epithelial", "epithelial", "epithelial", "tumor_microenv", "tumor_microenv", "stromal_rich", "stromal_rich", "stromal_rich", "stromal_rich", "lymphoid", "lymphoid", "lymphoid", "lymphoid", "lymphoid", "lymphoid", "brain_cortex", "brain_cortex", "brain_cortex", "brain_cortex", "brain_cortex", "columnar_muscle", "columnar_muscle", "columnar_muscle"),
  cell_type = c("epiteliale", "immunitario", "stromale e muscolo liscio", "epiteliale", "fibroblasti", "endoteliale", "epiteliale", "fibroblasti", "immunitario", "endoteliali", "linfociti B", "linfociti T", "mieloidi", "plasmacellule", "stromali", "astrociti", "microglia", "oligodendrociti", "vascolari", "neuroni", "cardiomiociti", "fibroblasti", "immunitari"),
  ell_um = c(32, 32, 8, 128, 32, 32, 128, 512, 16, NA, 512, 16, NA, NA, NA, 8, NA, 512, NA, 128, 64, 512, 8),
  s = c(5, 4, 6, 6, 4, 3.5, 6, 3, 2.5, 0, 4, 1, 0, 0, 0, 6, 0, 5, 0, 6, 3, 2, 5),
  nuc_median_area_um2 = c(13.9398, 11.2418, 11.2418, 32.512, 24.3465, 31.1635, 32.4371, 31.9876, 26.2942, 15.1405, 16.7894, 16.0399, 17.0143, 15.6651, 15.7401, 49.9341, 23.4277, 30.9366, 23.3902, 110.831, 19.7178, 19.9427, 17.9185),
  src = c("B-209", "B-210", "B-211", "B-212", "B-213", "B-214", "B-215", "B-216", "B-217", "B-218", "B-219", "B-220", "B-221", "B-222", "B-223", "B-224", "B-225", "B-226", "B-227", "B-228", "B-229", "B-230", "B-231"),
  stringsAsFactors = FALSE
)

.TP_SPATIAL_DEFAULT <- "uniform"   # S1.5c: cambia solo con la regola 3 pre-registrata

.TP_EPI_LABEL <- c(tumor_microenv = "tumorale", stromal_rich = "tumorale")   # DEC-014: istologia delle ROI (R2)

.tp_mut <- function() Sys.getenv("S15B_MUTANT", "none")

# Catalogo e composizione per tipo di un preset (S1.5b). Densita' del tipo puro rho_t = D * mbar / m_t (mbar = sum f_t m_t), cosi' la media
# armonica della miscela vale D; modello di taglia portato a rho_t con "preserve_area" (stessa funzione A_N(A_T)); braccio "recentred":
# intercetta spostata di (m_t / m_sim)^beta.
.tp_multi <- function(p, size_arm, composition) {
  ty <- .TP_TYPES[.TP_TYPES$preset == p$preset, , drop = FALSE]
  mut <- .tp_mut()
  fcol <- switch(composition, tissue = "fraction", follicle = "fraction_follicle", paracortex = "fraction_paracortex")
  f <- ty[[fcol]]
  if (anyNA(f)) stop("tissue_preset_spec(): composizione '", composition, "' non disponibile per il preset ", p$preset, call. = FALSE)
  if (mut == "M2") f <- f * 0.9                                                     # mutante: frazioni non rinormalizzate
  if (abs(sum(f) - 1) > 1e-5) stop("tissue_preset_spec(): le frazioni del preset non sommano a 1 (", format(sum(f), digits = 8), ")", call. = FALSE)
  f <- f / sum(f)                                                                     # arrotondamento a 6 cifre della tabella
  rho_g <- ty$rho_g; beta_g <- ty$beta_g; sigma_g <- ty$sigma_g
  if (mut == "M1" && p$preset == "brain_cortex") {                                  # mutante: modelli di neuroni e glia scambiati
    gN <- ty$size_group[ty$cell_type == "neuroni"]; iN <- ty$size_group == gN; iG <- !iN
    tmp <- c(rho_g[iN][1], beta_g[iN][1], sigma_g[iN][1])
    rho_g[iN] <- rho_g[iG][1]; beta_g[iN] <- beta_g[iG][1]; sigma_g[iN] <- sigma_g[iG][1]
    rho_g[iG] <- tmp[1]; beta_g[iG] <- tmp[2]; sigma_g[iG] <- tmp[3]
  }
  if (mut == "M6") { rho_g[] <- p$rho; beta_g[] <- p$beta; sigma_g[] <- p$sigma }   # mutante: gruppi di taglia ignorati
  mbar <- sum(f * ty$m_t)                                                             # composizione scelta: la miscela vale D
  d0 <- mbar / ty$m_t
  dens <- if (mut == "M3") p$density * d0 / sum(f * d0) else p$density * d0           # M3: media aritmetica invece di armonica
  conv <- if (mut == "M7") 1 else (dens / p$density_model)^((1 - beta_g) / 2)
  rho <- rho_g * conv
  if (size_arm == "recentred") rho <- rho * (if (mut == "M4") ty$m_sim / ty$m_t else ty$m_t / ty$m_sim)^(beta_g / 2)
  epi <- ifelse(ty$cell_type == "epiteliale" & p$preset %in% names(.TP_EPI_LABEL), .TP_EPI_LABEL[p$preset], NA_character_)
  if (mut == "M5") epi[!is.na(epi)] <- "normale"
  ct <- data.frame(cell_type = ty$cell_type, density = dens, min_dist_um = p$min_dist_um, nucleus_shape = "circle",
                   nucleus_radius = NA_real_, nucleus_to_eq_ratio = rho, nucleus_area_slope = beta_g, nucleus_area_sdlog = sigma_g,
                   nucleus_aspect = 1, nucleus_orient = "random", nucleus_nc_max = p$phi,
                   pale_fraction = 1 - p$density_evident / p$density, size_group = as.integer(ty$size_group),
                   epithelium_label = unname(epi), stringsAsFactors = FALSE)
  list(cell_types = ct, composition = data.frame(cell_type = ty$cell_type, fraction = f, stringsAsFactors = FALSE),
       type_sources = data.frame(cell_type = ty$cell_type, source = ty$src, stringsAsFactors = FALSE))
}

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
#' @param size_conversion `"preserve_area"` (predefinito) o `"none"` (solo con `types = "single"`).
#' @param types `"multi"` (predefinito da S1.5b: tipi cellulari) o `"single"` (un tipo per tessuto, preset di S1.5, invariato).
#' @param size_arm con `types = "multi"`: `"recentred"` (predefinito, NULL; intercetta per tipo ricentrata sul territorio simulato
#'   tipico, vincitore della regola di S1.5b) o `"group"` (modello del gruppo di taglia com'e' stimato).
#' @param composition `"tissue"`; per il linfonodo anche `"follicle"` e `"paracortex"` (DEC-015).
#' @param spatial con `types = "multi"`: `"uniform"` (semina cieca al tipo di S1.5b) o `"compartments"` (S1.5c: tabella
#'   `compartments` del preset, calibrata su tutte le ROI). NULL = predefinito del pacchetto (`.TP_SPATIAL_DEFAULT`).
#' @param min_dist_mode `"tissue"` (d* del tessuto per tutti i tipi, DEC-009) o `"nucleus"` (S1.5c: d_t = kappa * 2 sqrt(a_t / pi),
#'   kappa tale che sum f_t d_t = d*; solo con `spatial = "compartments"`).
#' @return lista con `cell_types` (catalogo: una riga, o una per tipo con `size_group` ed `epithelium_label`), `composition` (con
#'   `types = "multi"`: tipo, frazione), `type_sources` (fonti B-nnn per tipo), `c8_interval` (µm), `archetype`, `tissue`, `sources`
#'   (data.frame parametro, valore, fonte), `size_conversion`.
#' @export
tissue_preset_spec <- function(preset, size_conversion = c("preserve_area", "none"), types = c("multi", "single"),
                               size_arm = NULL, composition = c("tissue", "follicle", "paracortex"), spatial = NULL,
                               min_dist_mode = c("tissue", "nucleus")) {
  size_conversion <- match.arg(size_conversion); types <- match.arg(types); composition <- match.arg(composition)
  min_dist_mode <- match.arg(min_dist_mode)
  if (is.null(spatial)) spatial <- .TP_SPATIAL_DEFAULT
  spatial <- match.arg(spatial, c("uniform", "compartments"))
  if (types == "single" && (spatial == "compartments" || min_dist_mode == "nucleus"))
    stop("tissue_preset_spec(): compartimenti e distanza per nucleo richiedono types = \"multi\" (S1.5c).", call. = FALSE)
  if (min_dist_mode == "nucleus" && spatial != "compartments")
    stop("tissue_preset_spec(): min_dist_mode = \"nucleus\" richiede spatial = \"compartments\" (S1.5c).", call. = FALSE)
  if (types == "multi") {
    if (is.null(size_arm)) size_arm <- "recentred"   # braccio vincente della regola pre-registrata di S1.5b (results/S1.5b/S1.5b_arm_rule.csv)
    if (!size_arm %in% c("group", "recentred"))
      stop("tissue_preset_spec(): con types = \"multi\" serve size_arm = \"group\" o \"recentred\" (S1.5b).", call. = FALSE)
    if (size_conversion != "preserve_area") stop("tissue_preset_spec(): con types = \"multi\" la conversione e' sempre \"preserve_area\".", call. = FALSE)
  }
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
  out <- list(cell_types = ct, c8_interval = c(p$c8_lo, p$c8_hi), archetype = p$archetype, tissue = p$tissue,
              sources = src, size_conversion = size_conversion, preset = preset, types = types, size_arm = NA_character_,
              composition = NULL, composition_name = NA_character_, spatial = if (types == "multi") spatial else "uniform",
              min_dist_mode = min_dist_mode, compartments = NULL)
  if (types == "multi") {
    m <- .tp_multi(p, size_arm, composition)
    out$cell_types <- m$cell_types; out$composition <- m$composition; out$type_sources <- m$type_sources
    out$size_arm <- size_arm; out$composition_name <- composition
    cp <- .TP_COMPART[.TP_COMPART$preset == preset, , drop = FALSE]
    cp <- cp[match(out$cell_types$cell_type, cp$cell_type), , drop = FALSE]
    out$compartments <- data.frame(cell_type = cp$cell_type, ell_um = cp$ell_um, s = cp$s, source = cp$src, stringsAsFactors = FALSE)
    if (min_dist_mode == "nucleus") {
      f <- out$composition$fraction[match(out$cell_types$cell_type, out$composition$cell_type)]
      raw <- 2 * sqrt(cp$nuc_median_area_um2 / pi); kappa <- p$min_dist_um / sum(f * raw)
      out$cell_types$min_dist_um <- kappa * raw; out$min_dist_kappa <- kappa
    }
  }
  out
}
