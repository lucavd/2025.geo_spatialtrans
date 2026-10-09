# results/S1.3/review/claims/02_diagnostics.R — verifiche mirate di affermazioni del report (dopo 01_recompute.R)
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
OD <- "results/S1.3/review/claims"; TOL <- 1e-9
c5 <- read.csv(file.path(OD, "rc_c5_from_cells.csv")); rr <- read.csv(file.path(OD, "rc_roi_real.csv"))
cat("== C-5 dalle celle (geos) vs vettore memorizzato\n")
g <- rr[rr$engine == "geos", ]; m <- merge(c5, g, by = c("archetype", "roi_id", "method"))
cat("interior_identical concordi:", sum(m$my_interior_identical == (m$interior_identical == 1)), "/40; n_interior_diff concordi:", sum(m$my_n_interior_diff == m$n_interior_diff),
    "; max|rel_area_interior mia-mem|:", max(abs(m$my_rel_area_interior - m$rel_area_interior)), "; max|rel_area_mono mia-mem|:", max(abs(m$my_rel_area_mono - m$rel_area_mono)),
    "; max|rel_summary mia-mem|:", max(abs(m$my_rel_summary_vs_R3 - m$rel_summary)), "; max mio riassunto vs memorizzato:", max(m$my_summary_vs_stored), "\n")
cat("ROI con flag d'interna diversi su celle SENZA frammenti:", sum(c5$n_interior_diff_mono > 0), "; celle:", sum(c5$n_interior_diff_mono), "; su celle con frammenti:", sum(c5$my_n_interior_diff - c5$n_interior_diff_mono), "\n")
cat("ROI con aree diverse su interne senza frammenti (>1e-9):", sum(c5$my_rel_area_interior_mono > TOL), "; max", max(c5$my_rel_area_interior_mono), "\n")
print(c5[c5$my_rel_area_mono > TOL | c5$n_interior_diff_mono > 0, c("archetype", "roi_id", "method", "my_rel_area_mono", "my_rel_area_interior_mono", "n_interior_diff_mono", "nsides_eq_int_mono", "max_abs_decc_int_mono")], row.names = FALSE)
cat("nsides identici sulle interne senza frammenti: min frazione", min(c5$nsides_eq_int_mono), "; max |Δecc|", max(c5$max_abs_decc_int_mono), "\n")
cat("ROI con rel_summary > 1e-9 ma C-5 testo PASS:\n"); g$c5_text <- as.logical(g$c5_text); print(g[g$c5_text & g$rel_summary > TOL, c("archetype", "roi_id", "method", "rel_summary")], row.names = FALSE)
d <- rr[rr$engine == "deldir", ]; cat("deldir ROI con rel_area_mono > 1e-9:\n"); print(d[d$rel_area_mono > TOL, c("archetype", "roi_id", "method", "rel_area_mono", "n_snapped")], row.names = FALSE)
cat("geos ROI con rel_area_mono > 1e-9:\n"); print(g[g$rel_area_mono > TOL, c("archetype", "roi_id", "method", "rel_area_mono", "n_snapped")], row.names = FALSE)
arb <- read.csv(file.path(OD, "rc_arbitration.csv"))
cat("arbitro: discordanze per ROI reale:\n"); print(table(paste(arb$set, arb$archetype, arb$roi_id, arb$model)))
cat("arbitro, min_edge > 2e-5:\n"); print(arb[arb$min_edge > 2e-5, c("set", "archetype", "roi_id", "model", "cell", "a_geos", "a_deldir", "a_hp", "ns_geos", "ns_deldir", "ns_hp", "min_edge", "geos_ok", "deldir_ok")], row.names = FALSE)
cat("quantili min_edge:", signif(quantile(arb$min_edge, c(0, .5, .9, .95, 1)), 3), "\n")
c6 <- read.csv(file.path(OD, "rc_c6.csv")); nf <- c6[c6$n_fragments == 0 & !c6$ident, ]
cat("== C-6: repliche senza frammenti NON identiche (", nrow(nf), ")\n"); print(nf, row.names = FALSE)
cat("di queste, con replica 1 (arbitrata):", sum(nf$rep == 1), "\n")
cat("con frammenti: worst col tab\n"); print(table(c6$worst[c6$n_fragments > 0 & !c6$ident]))
cat("== A5 r5: area della regione e C-1\n")
rs <- read.csv("results/R3/R3_roi_summary.csv"); a <- rs$area_poly_um2[rs$archetype == "A5" & rs$roi_id == "r5" & rs$method == "cellpose_rgb"]
ns <- read.csv("results/S1.3/S1.3_null_summary.csv"); x <- ns[ns$archetype == "A5" & ns$roi_id == "r5" & ns$model == "RSA" & ns$rep == 15, ]
cv <- read.csv("results/S1.3/S1.3_c2_coverage.csv"); y <- cv[cv$archetype == "A5" & cv$roi_id == "r5", ]
cat("area_poly", a, "; c1 rel", x$c1, "-> abs", x$c1 * a, "; c2_union rel", y$c2_union, "-> abs", y$c2_union * a, "; c2_overlap", x$c2_overlap, "c2_symdiff", x$c2_symdiff, "->abs", x$c2_symdiff * a, "\n")
cat("== nulli geos con C-2 FAIL: elenco\n"); f <- ns[ns$c2_overlap > TOL | ns$c2_symdiff > TOL, c("archetype", "roi_id", "model", "rep", "c2_overlap", "c2_symdiff", "n_multipart")]; print(f, row.names = FALSE)
cat("== checkB: B3 G2 A5 range e numero ROI entro 10 %\n"); b <- read.csv(file.path(OD, "rc_checkB_roi.csv"))
for (md in c("RSA")) for (A in paste0("A", 1:5)) { z <- b[b$archetype == A & b$model == md, ]; cat(A, md, "B2", sprintf("%+.3f..%+.3f", min(z$B2), max(z$B2)), "B3", sprintf("%+.3f..%+.3f", min(z$B3), max(z$B3)), "n|B3|<=.1:", sum(abs(z$B3) <= .1), "\n") }
ks <- read.csv(file.path(OD, "rc_ks_nsides.csv")); print(ks[ks$archetype %in% c("A4", "A6"), ], row.names = FALSE)
t <- read.csv("results/S1.3/S1.3_test_deldir_none.csv"); print(t[t$status == "FAIL", c("check", "input", "metric", "value")], row.names = FALSE)
t4 <- read.csv("results/S1.3/S1.3_test_geos_M4.csv"); cat("M4 FAIL metrics:\n"); print(table(t4$metric[t4$status == "FAIL"]))
