# R/testing/test_S2.8_d.R - S2.8, gruppo d della batteria (letto da test_S2.8.R): test aggiunti dopo la mappa del documento di progetto
# (results/S2.8/S2.8_design_map.csv). Valori della tabella dei parametri (design sezione 9) alla precisione scritta nel documento,
# cuore con la mediana delle altre sezioni (D2-5), tumore e stroma sulla stessa sezione, fondo dei due motori da capture_params(),
# argomenti della funzione a pezzi senza test (ceiling_q, assemble, keep_profiles, archetype, verbose).
TAB <- list(   # valori come scritti nel documento (stringhe: la tolleranza e' mezza unita' dell'ultima cifra scritta)
  r = c("0.75", "0.75", "0.82", "0.55", "1.48", "2.61"), sigma_e = c("6.23", "4.72", "4.72", "4.49", "4.99", "4.86"),
  D = c("24.75", "8.94", "8.94", "22.04", "6.86", "15.49"), w = c("0.087", "0.040", "0.040", "0.109", "0.112", "0.098"),
  ell = c("32.85", "24.73", "24.73", "41.50", "33.07", "32.96"), field_s2 = c("0.5", "0.5", "0.354", "0.297", "0.0743", "0.841"),
  field_L = c("0.6375", "0.5", "0.85", "1.275", "1.275", "0.885"))
CAP <- list(b = c("0", "0.174", "0.174", "0", "0.147", "0"), parity_row = c("1.417", "1.578", "1.578", "0.996", "1.201", "0.947"),
  parity_col = c("1.750", "1.724", "1.724", "1.105", "1.231", "0.967"), sd_alpha = c("0.193", "0.221", "0.221", "0.211", "0.242", "0.183"),
  sd_beta = c("0.185", "0.191", "0.191", "0.251", "0.180", "0.188"), rho_alpha = c("-0.305", "-0.313", "-0.313", "-0.390", "-0.310", "-0.174"),
  rho_beta = c("-0.275", "-0.212", "-0.212", "-0.461", "-0.296", "-0.208"), w2 = c("0.0362", "0.0232", "0.0232", "0.0363", "0.0412", "0.0774"),
  ell2 = c("1826", "429", "429", "1575", "305", "113"))
tol_of <- function(x) { d <- if (grepl("\\.", x)) nchar(sub(".*\\.", "", x)) else 0L; 0.5 * 10^(-d) + 1e-12 }
AR <- paste0("A", 1:6)
chk_tab <- function(tab, getter, label) for (nm in names(tab)) {
  v <- vapply(AR, function(A) getter(A)[[nm]], numeric(1)); e <- as.numeric(tab[[nm]]); tl <- vapply(tab[[nm]], tol_of, numeric(1))
  ok(sprintf("d1 %s `%s` = tabella della sezione 9 alla precisione scritta (scarto massimo %.2g)", label, nm, max(abs(v - e))), all(abs(v - e) <= tl))
}
chk_tab(TAB, transcript_params, "transcript_params()"); chk_tab(CAP, capture_params, "capture_params()")
ok("d1 sigma_in = 0.7 um in tutti i tessuti (DEC-038)", all(vapply(AR, function(A) transcript_params(A)$sigma_in, numeric(1)) == 0.7))
oth <- c("A1", "A2", "A4", "A5"); tp <- lapply(AR, transcript_params); names(tp) <- AR
ok("d2 cuore: sigma_e, D, w, ell = mediane delle altre sezioni (D2-5, DEC-037; entro 1e-3 relativo)",
   all(vapply(c("sigma_e", "D", "w", "ell"), function(k) abs(tp$A6[[k]] / stats::median(vapply(oth, function(A) tp[[A]][[k]], numeric(1))) - 1) < 1e-3, logical(1))))
ok("d3 tumore e stroma: stesso spostamento (sigma_e, D, w, ell) e stessa cattura (sezione A23); r e campo diversi",
   identical(tp$A2[c("sigma_e", "D", "w", "ell")], tp$A3[c("sigma_e", "D", "w", "ell")]) &&
   identical(capture_params("A2")[setdiff(names(capture_params("A2")), "archetype")], capture_params("A3")[setdiff(names(capture_params("A3")), "archetype")]) &&
   tp$A2$r != tp$A3$r && tp$A2$field_L != tp$A3$field_L)
# fondo: i due motori usano capture_params()$b, non transcript_params()$b (M-106); linfonodo con b = 0 -> si prova il tumore
L2 <- mk_layer("tumor_microenv", 100L, 202L)
op <- simulate_visium_hd_pieces(L2, seed = 3L, piece_um = Inf, frame_open = FALSE); os5 <- simulate_visium_hd(L2, seed = 3L, frame_open = FALSE)
ok(sprintf("d4 fondo dei due motori = capture_params()$b (%.4f), non quello di R5 (%.4f)", capture_params("A2")$b, transcript_params("A2")$b),
   identical(op$truth$layers$background$b, capture_params("A2")$b) && identical(os5$truth$layers$background$b, capture_params("A2")$b))
# argomenti della funzione a pezzi senza test (mappa, S2.8_exported_args.csv)
base <- function(...) simulate_visium_hd_pieces(L4, seed = 7L, piece_um = 112, frame_open = TRUE, ...)
o0 <- base(); oq <- base(ceiling_q = 1e-3)
ok(sprintf("d5 ceiling_q: tetto = stripe_ceiling(q) (%.4f con 1e-3 contro %.4f con 1e-5) e conteggi diversi", oq$metadata$overgeneration, o0$metadata$overgeneration),
   identical(oq$metadata$overgeneration, stripe_ceiling(capture_params("A4"), 1e-3)) && oq$metadata$overgeneration < o0$metadata$overgeneration &&
   !identical(oq$counts_2um, o0$counts_2um))
ol <- base(assemble = "light")
ok("d6 assemble = \"light\": UMI per quadratino = colonne di counts_2um, counts_8um e celle identici al completo",
   identical(as.numeric(ol$umi_2um), as.numeric(Matrix::colSums(o0$counts_2um))) && identical(ol$counts_8um, o0$counts_8um) &&
   identical(ol$truth$cells, o0$truth$cells) && is.null(ol$counts_2um))
ok_ <- base(keep_profiles = FALSE)
ok("d7 keep_profiles = FALSE: nessuna matrice dei profili, conteggi identici", is.null(ok_$truth$cell_counts) && !is.null(o0$truth$cell_counts) &&
   identical(ok_$counts_2um, o0$counts_2um))
oa <- base(archetype = "A4"); ov <- suppressMessages(base(verbose = TRUE))
ok("d8 archetype esplicito = archetipo dal preset; verbose = TRUE non cambia nulla", identical(oa$counts_2um, o0$counts_2um) && identical(ov$counts_2um, o0$counts_2um))
ok("d9 archetype sbagliato: errore chiaro", grepl("archetipo sconosciuto", err_of(base(archetype = "A9"))))
