"""R/testing/mutants_S1.6.py - mutanti N1-N6 di S1.6 (docs/steps/S1.6.md, G-2): ognuno deve far fallire almeno un'asserzione di test_S1.6.R."""
P2 = "R/04b2_seed_centroids.R"; P5 = "R/04b5_cell_layer_validation.R"; PC = "R/04b_cell_placement.R"
WM = "    wmax <- if (b_used == \"exact\") max(rho) * hmax * (1 + 1e-12) else 1.1 * max(wg)"
MUTANTS = {
    "N1": [(P2, "quotas = getOption(\"geo.spatialtrans.quotas\", FALSE),", "quotas = getOption(\"geo.spatialtrans.quotas\", TRUE),")],   # quote predefinite
    "N2": [(P2, WM, "    wmax <- 1.1 * max(wg)")],                                                                                # limite di oggi
    "N3": [(P2, WM, "    wmax <- if (b_used == \"exact\") 0.5 * max(rho) * hmax else 1.1 * max(wg)")],                         # limite troppo basso
    "N4": [(P5, "    blind <- identical(type_order, \"blind\") && \"target_density_weighted\" %in% names(rd)", "    blind <- FALSE")],   # C7 di prima
    "N5": [(PC, "  t_order <- if (!is.null(type_order)) type_order else if", "  t_order <- if (FALSE) type_order else if")],    # type_order ignorato
    "N6": [("R/testing/test_S1.1.R", "write.csv(TT, file.path(OUTD, \"S1.1_test_results.csv\")", "write.csv(TT, file.path(RES, \"S1.1_test_results.csv\")")],
}
