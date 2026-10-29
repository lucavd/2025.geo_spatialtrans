"""R/testing/mutants_S1.6.py - mutanti di S1.6 (docs/steps/S1.6.md, G-2): ognuno deve far fallire almeno un'asserzione di test_S1.6.R.
N1-N6 pre-registrati; N7-N16 sono i mutanti R1-R10 del revisore del codice (G-5, /mnt/micron/geo_spatialtrans/S1.6/review_code/mutants_review.py),
aggiunti dopo la revisione con i test che li rilevano (T1c, T2e, T2f, T3c, T3d, T4d-T4f, T7)."""
P2 = "R/04b2_seed_centroids.R"; P5 = "R/04b5_cell_layer_validation.R"; PC = "R/04b_cell_placement.R"
WM = "    wmax <- if (b_used == \"exact\") max(rho) * hmax * (1 + 1e-12) else 1.1 * max(wg)"
MUTANTS = {
    "N1": [(P2, "quotas = getOption(\"geo.spatialtrans.quotas\", FALSE),", "quotas = getOption(\"geo.spatialtrans.quotas\", TRUE),")],   # quote predefinite
    "N2": [(P2, WM, "    wmax <- 1.1 * max(wg)")],                                                                                # limite di oggi
    "N3": [(P2, WM, "    wmax <- if (b_used == \"exact\") 0.5 * max(rho) * hmax else 1.1 * max(wg)")],                         # limite troppo basso
    "N4": [(P5, "    blind <- identical(type_order, \"blind\") && \"target_density_weighted\" %in% names(rd)", "    blind <- FALSE")],   # C7 di prima
    "N5": [(PC, "  t_order <- if (!is.null(type_order)) type_order else if", "  t_order <- if (FALSE) type_order else if")],    # type_order ignorato
    "N6": [("R/testing/test_S1.1.R", "write.csv(TT, file.path(TOUT, \"S1.1_test_results.csv\")", "write.csv(TT, file.path(RES, \"S1.1_test_results.csv\")")],
    # dal revisore (R1-R10)
    "N7": [(P5, "    blind <- identical(type_order, \"blind\") && ", "    blind <- !identical(type_order, \"sequential\") && ")],
    "N8": [(P5, "    if (is.null(type_order)) type_order <- cl$metadata$seed_info$type_order\n", "    if (is.null(type_order)) type_order <- \"sequential\"\n")],
    "N9": [(PC, "type_order = if (sp_used == \"compartments\") \"compartments\" else t_order)", "type_order = t_order)")],
    "N10": [(PC, "quotas = quotas, thinning_bound = thinning_bound, verbose = FALSE)", "quotas = quotas, verbose = FALSE)")],
    "N11": [(PC, "heterogeneity = htab, quotas = quotas, thinning_bound", "heterogeneity = htab, thinning_bound")],
    "N12": [(P2, "      hg <- NULL; a <- .sc_solve_a(eta_g, f, rho)\n", "      hg <- NULL; a <- .sc_solve_a(eta_g, f, rho); hmax <- NA_real_\n")],
    "N13": [(P2, "max(rho) * hmax * (1 + 1e-12)", "max(type_df$density) * hmax * (1 + 1e-12)")],
    "N14": [(PC, "types = getOption(\"geo.spatialtrans.preset_types\", \"multi\"), size_arm = NULL,", "types = \"multi\", size_arm = NULL,")],
    "N15": [(P2, "else if (is.na(cap$s)) 1 else cap$H   #", "else 1   #")],
    "N16": [("R/testing/test_S1.2.R", "write.csv(res, file.path(TOUT, \"S1.2_test_results.csv\"), row.names = FALSE)",
             "write.csv(res, paste0(RES, \"/S1.2_test_results.csv\"), row.names = FALSE)")],
}
