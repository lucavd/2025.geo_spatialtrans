"""R/testing/mutants_S1.5.py - mutanti M1-M9 di S1.5 (stessi comportamenti degli interruttori S15_MUTANT tolti dal codice in S1.6, BL-155)."""
P5 = "R/04b5_cell_layer_validation.R"; PC = "R/04b_cell_placement.R"; P0 = "R/04b0_tissue_presets.R"; P6 = "R/04b6_estimate_pixel_size.R"
A21 = "    a21 <- { z <- geos::geos_area(geos::geos_intersection_prec(pT_loc[pr[, 2]], pT_loc[pr[, 1]], g)); n_eval <- n_eval + length(z); z }\n"
MUTANTS = {
    # M1: partizione misurata con st_union (BL-070)
    "M1": [(P5, "  nb2 <- nb_sum + nb_pair + nb_out\n",
            "  U <- vapply(rd$region_id, function(r) as.numeric(sf::st_area(sf::st_union(Tt[cd$region_id == r]))), 0)\n"
            "  nb_sum <- sum(abs(U - aR) / aR > 0.01); nb_pair <- 0L; nb_out <- 0L\n  nb2 <- nb_sum + nb_pair + nb_out\n")],
    # M2: una sola direzione dell'intersezione delle coppie
    "M2": [(P5, A21, "    a21 <- a12\n")],
    # M3: NA trattati come 0 nella somma e nelle coppie
    "M3": [(P5, "  rel_sum <- (abs(as.numeric(sT) + gap - aR) / aR)[has_cells]\n",
            "  rel_sum <- (abs(as.numeric(sT) + gap - aR) / aR)[has_cells]\n  rel_sum[is.na(rel_sum)] <- 0\n"),
           (P5, A21, A21 + "    a12[is.na(a12)] <- 0; a21[is.na(a21)] <- 0\n")],
    # M4: stesso seme per la semina e per la taglia
    "M4": [(PC, "  s <- sample.int(.Machine$integer.max, 3L)\n", "  s <- sample.int(.Machine$integer.max, 3L)\n  s[2] <- s[1]\n")],
    # M5: nuclei permutati rispetto alle cellule
    "M5": [(PC, "  nuc <- dg$cell_nuclei\n", "  nuc <- dg$cell_nuclei\n  if (length(nuc) > 1) nuc <- nuc[c(2:length(nuc), 1)]\n")],
    # M6: phi di linfonodo e cuore scambiati
    "M6": [(P0, "  tb <- .TP_TABLE\n", "  tb <- .TP_TABLE\n  tb$phi[c(4, 6)] <- tb$phi[c(6, 4)]\n")],
    # M7: scala ignorata (1 um/px)
    "M7": [(PC, "  px_used <- px\n", "  px_used <- 1\n")],
    # M8: tag di risoluzione del TIFF (96 dpi) invece dei metadati (BL-016)
    "M8": [(P6, "  j <- jsonlite::fromJSON(f)\n",
            "  return(list(pixel_size_um = 25400 / 96, method = \"metadata\", n_features = NA_integer_, ci = c(NA_real_, NA_real_),\n"
            "              sigma_px = NA_real_, refused = FALSE, reason = \"\", detail = list()))\n  j <- jsonlite::fromJSON(f)\n")],
    # M9: rifiuti della stima della scala disattivati
    "M9": [(P6, "  if (s_star < 1) return(refuse(", "  if (FALSE && s_star < 1) return(refuse("),
           (P6, "  if (kstar == length(sig)) return(refuse(", "  if (FALSE && kstar == length(sig)) return(refuse("),
           (P6, "  if (N < 30) return(refuse(", "  if (FALSE && N < 30) return(refuse("),
           (P6, "  if (size_check && is.finite(r_exp) && (r_obs", "  if (FALSE && size_check && is.finite(r_exp) && (r_obs")],
}
