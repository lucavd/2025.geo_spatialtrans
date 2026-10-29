"""R/testing/mutants_S1.5b.py - mutanti M1-M8 di S1.5b (stessi comportamenti degli interruttori S15B_MUTANT tolti dal codice in S1.6, BL-155)."""
P0 = "R/04b0_tissue_presets.R"; PC = "R/04b_cell_placement.R"
RG = "  rho_g <- ty$rho_g; beta_g <- ty$beta_g; sigma_g <- ty$sigma_g\n"
EPI = "  epi <- ifelse(ty$cell_type == \"epiteliale\" & p$preset %in% names(.TP_EPI_LABEL), .TP_EPI_LABEL[p$preset], NA_character_)\n"
MUTANTS = {
    # M1: modelli di taglia di neuroni e glia scambiati nel cervello
    "M1": [(P0, RG, RG + "  if (p$preset == \"brain_cortex\") {\n    gN <- ty$size_group[ty$cell_type == \"neuroni\"]; iN <- ty$size_group == gN; iG <- !iN\n"
            "    tmp <- c(rho_g[iN][1], beta_g[iN][1], sigma_g[iN][1])\n    rho_g[iN] <- rho_g[iG][1]; beta_g[iN] <- beta_g[iG][1]; sigma_g[iN] <- sigma_g[iG][1]\n"
            "    rho_g[iG] <- tmp[1]; beta_g[iG] <- tmp[2]; sigma_g[iG] <- tmp[3]\n  }\n")],
    # M2: frazioni non rinormalizzate
    "M2": [(P0, "  if (abs(sum(f) - 1) > 1e-5)", "  f <- f * 0.9\n  if (abs(sum(f) - 1) > 1e-5)")],
    # M3: media aritmetica invece di armonica
    "M3": [(P0, "  dens <- p$density * d0\n", "  dens <- p$density * d0 / sum(f * d0)\n")],
    # M4: ricentratura invertita
    "M4": [(P0, "(ty$m_t / ty$m_sim)^(beta_g / 2)", "(ty$m_sim / ty$m_t)^(beta_g / 2)")],
    # M5: etichetta dell'epitelio sempre normale
    "M5": [(P0, EPI, EPI + "  epi[!is.na(epi)] <- \"normale\"\n")],
    # M6: gruppi di taglia ignorati
    "M6": [(P0, RG, RG + "  rho_g[] <- p$rho; beta_g[] <- p$beta; sigma_g[] <- p$sigma\n")],
    # M7: conversione del modello di taglia ignorata
    "M7": [(P0, "  conv <- (dens / p$density_model)^((1 - beta_g) / 2)\n", "  conv <- 1\n")],
    # M8: semina in sequenza invece che cieca per i preset per tipo
    "M8": [(PC, "else if (!is.null(spec) && identical(spec$types, \"multi\")) \"blind\" else \"sequential\"", "else \"sequential\"")],
}
