# test_arm_p_o.R — arm P (simulated μ^P, φ by rank, simulated method reference) and arm O (oracle μ)
# on the colon pilot reference (R/arm_p.R, R/simulate_spots.R; ADEMP §2.1; DECISIONS 2026-11-04).
# Pilot data only (<TESI_DATA>/derived/reference/colon_pilot.rds: pilot held-out donor drawn inside
# training), plus the main P calibration (colon_p_calibration.rds); no main held-out cell is read.
#
# Pre-registered checks (written 2026-11-05 before the first run; PASS/FAIL printed, stop at the end
# if any FAIL). One replicate, seed 42, 500 spots, 10 cells per spot, 4,900 UMI per cell (the level
# where differences between μ are least hidden by sampling noise).
#   P1. μ^P: genes = reference genes (named P00001...), the 18 reference types, no NA, every column
#       sums to 1 (|1 - sum| < 1e-12); markers per type = the calibration counts.
#   P2. Reproducible: same seed -> identical μ^P; new seed -> different μ^P.
#   P3. No two P types have the same profile (all types have markers; Pearson of every pair < 1 - 1e-6).
#   P4. φ^P: the same set of values as the real φ (sorted vectors identical); the rank of φ^P over the
#       P genes reproduces the rank of φ over the real genes (Spearman of (base mean, φ^P) equal to
#       Spearman of (pooled real mean, φ) within 0.01: only ties can differ).
#   P5. P reference: same number of cells, types and donors as the real reference; column totals /
#       L: mean within 1% of 1; per type, Pearson of the mean simulated profile with μ^P ≥ 0.99.
#   O1. μ_O comes from one donor, the pilot held-out donor, none of whose cells is in the reference;
#       columns in (0.9, 1]; second route: μ_O recomputed from the saved held-out counts (row sums
#       per type, rescaled) is proportional to mu_held, max abs difference after normalising
#       both columns to 1 < 1e-12.
#   O2. Spots of E, C−, O, P and the anchor: 500 each, same compositions (one cell table), after
#       thinning no spot above the target; spots of P have the P gene names.
# Counter-proofs (expected results, reported; a result against the expectation is written, not hidden):
#   K1. Mean over spots of the Pearson correlation (log1p CP10k, kept genes) between spot s of an
#       arm and spot s of the anchor: O > E > C−.
#   K2. NNLS on the P spots: with the P reference, RMSE below the null baseline; with the real
#       reference (P genes relabelled as the real genes), RMSE at or above the null baseline, like
#       C− (the reason for Luca's decision of 2026-11-04).
# Reported: seconds for μ^P + φ^P + P reference + P spots (input to the cost check, ADEMP §6).
#
# Result of the first run (2026-11-05; expectations above unchanged): P1–P5, O1, O2 PASS; K1 as
# expected (O 0.877 > E 0.851 > C− 0.002); K2 AGAINST: NNLS RMSE 0.205 with the P reference, 0.200 with
# the real one, null 0.091. Part of the expectation was mis-specified: in the colon pilot NNLS does not
# beat the null baseline at 10 cells even on E (0.115 vs 0.090, colon_pilot_rmse_methods.csv), so
# "below the null" was not a fair bar. The finding that remains: with its own reference P scores
# like C− (0.234 in the pilot), not like E. Diagnosis (scratch, same seed): P markers sit at random
# on the gene-mean distribution (median percentile 0.47, vs 0.93 for the real markers), so a spot type
# puts 0.1–40% (median 1.7%) of its counts in its markers vs 2.5–60% (median 11.2%) in E, and P types
# differ little (median total-variation distance 0.076 vs 0.325 in E). Not a code error: markers on
# random genes are the definition of P (DECISIONS 2026-11-04, Luca; as in Splatter). P is kept as
# pre-registered (PROTOCOLLO §3.1) and its weak separability is reported in Aim 1. The real-reference
# part (0.200 ≈ C−) confirms Luca's reason for the simulated P reference. Pilot numbers, not results.
#
# Run from thesis/:  Rscript tests/test_arm_p_o.R

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
source("R/simulate_spots.R"); source("R/arm_p.R"); source("R/pilot_metrics.R"); source("R/deconvolve.R")
res <- character(0)
ok <- function(cond, msg) { r <- if (isTRUE(cond)) "PASS" else "FAIL"; res <<- c(res, r); cat(r, msg, "\n") }

R <- readRDS(file.path(data_dir, "derived", "reference", "colon_pilot.rds"))
calib <- readRDS(file.path(data_dir, "derived", "p_calibration", "colon_p_calibration.rds"))
G <- nrow(R$mu_E); types <- R$types_ref; ts <- R$types_spot
n_spots <- 500; n_cells <- 10; target <- n_cells * 4900
set.seed(42)

# --- arm P ------------------------------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
s0 <- .Random.seed
P <- simulate_mu_p(calib, G, types)
ok(identical(dim(P$mu), c(G, length(types))) && all(grepl("^P\\d{5}$", rownames(P$mu))) &&
   identical(colnames(P$mu), types) && !anyNA(P$mu) && all(abs(colSums(P$mu) - 1) < 1e-12) &&
   all(lengths(P$markers)[types] == lengths(calib$markers)[types]), "P1 μ^P dimensions, columns sum to 1, markers per type")
s1 <- .Random.seed
assign(".Random.seed", s0, envir = .GlobalEnv); P_again <- simulate_mu_p(calib, G, types)
P_new <- simulate_mu_p(calib, G, types)
assign(".Random.seed", s1, envir = .GlobalEnv)
ok(identical(P_again$mu, P$mu) && !identical(P_new$mu, P$mu), "P2 same seed same μ^P, new seed new μ^P")
cc <- cor(P$mu); diag(cc) <- 0
ok(max(cc) < 1 - 1e-6, sprintf("P3 no two P types identical (max pairwise Pearson %.4f)", max(cc)))
mean_real <- pooled_mean(R$mu_E, R$L_by_type)
phi_p <- phi_by_rank(R$phi, mean_real, P$base)
sp_p <- cor(P$base, phi_p, method = "spearman"); sp_r <- cor(mean_real, R$phi, method = "spearman")
ok(identical(sort(unname(phi_p)), sort(unname(R$phi))) && abs(sp_p - sp_r) < 0.01,
   sprintf("P4 φ^P = real φ values, Spearman(mean, φ) P %.4f vs real %.4f", sp_p, sp_r))
L_ref <- Matrix::colSums(R$ref)
RP <- simulate_reference_p(P$mu, phi_p, R$ref_meta, L_ref, colnames(R$ref))
ratio <- Matrix::colSums(RP$ref) / L_ref
prof_cor <- vapply(types, function(k) cor(RP$profiles[, k], P$mu[, k]), numeric(1))
ok(ncol(RP$ref) == ncol(R$ref) && identical(RP$ref_meta, R$ref_meta) && abs(mean(ratio) - 1) < 0.01 &&
   all(prof_cor >= 0.99), sprintf("P5 P reference: %d cells, totals/L mean %.4f, min profile Pearson %.4f (%s)",
   ncol(RP$ref), mean(ratio), min(prof_cor), names(which.min(prof_cor))))
t_p <- proc.time()[["elapsed"]] - t0

# --- arm O ------------------------------------------------------------------------------------------
mu_O <- mu_oracle(R$mu_held, ts)
hd <- unique(R$held_meta$donor_id)
re_mu <- sapply(ts, function(k) Matrix::rowSums(R$held[, R$held_meta$cell_type == k, drop = FALSE]))
d_route <- max(abs(sweep(re_mu, 2, colSums(re_mu), "/") - sweep(mu_O, 2, colSums(mu_O), "/")))
ok(length(hd) == 1 && hd == R$donors$heldout && !any(R$ref_meta$donor_id == hd) &&
   all(colSums(mu_O) > 0.9 & colSums(mu_O) <= 1) && d_route < 1e-12,
   sprintf("O1 μ_O from the held-out donor %s only; second route max diff %.1e", hd, d_route))

# --- one replicate ----------------------------------------------------------------------------------
n_train <- lengths(R$L_by_type[ts]); alpha <- length(n_train) * n_train / sum(n_train)
comp <- draw_compositions(n_spots, n_cells, ts, alpha)
cells <- draw_cells(comp, R$L_by_type[ts])
truth <- matrix(0, n_spots, length(types), dimnames = list(rownames(comp), types)); truth[, ts] <- comp / n_cells
held_ok <- R$held_meta$cell_type %in% ts
anc <- anchor_spots(R$held[, held_ok], R$held_meta$cell_type[held_ok], comp)
mu_E <- R$mu_E[, ts, drop = FALSE]
data <- list(E = simulate_arm(mu_E, R$phi, cells, rownames(comp)),
             Cminus = simulate_arm(permute_mu(mu_E), R$phi, cells, rownames(comp)),
             O = simulate_arm(mu_O, R$phi, cells, rownames(comp)),
             anchor = anc$counts)
t1 <- proc.time()[["elapsed"]]
data$P <- simulate_arm(P$mu[, ts, drop = FALSE], phi_p, cells, rownames(comp))
t_p <- t_p + proc.time()[["elapsed"]] - t1
data <- lapply(data, thin_to_depth, target = target)
ok(all(vapply(data, ncol, 1) == n_spots) && all(vapply(data, function(x) max(Matrix::colSums(x)), 1) <= target) &&
   identical(rownames(data$P), rownames(P$mu)) && identical(rownames(data$O), rownames(R$mu_E)),
   "O2 500 spots per dataset, same cells, none above the target, P spots on P genes")

# --- counter-proofs ---------------------------------------------------------------------------------
lcp <- function(x) { x <- as.matrix(x); log1p(sweep(x, 2, colSums(x), "/") * 1e4) }
la <- lcp(data$anchor)
k1 <- vapply(c("O", "E", "Cminus"), function(d) { lx <- lcp(data[[d]])
  mean(vapply(seq_len(n_spots), function(s) cor(lx[, s], la[, s]), 1)) }, 1)
cat(sprintf("K1 mean per-spot Pearson with the anchor: O %.4f, E %.4f, C− %.4f -> expected O > E > C−: %s\n",
            k1["O"], k1["E"], k1["Cminus"], if (k1["O"] > k1["E"] && k1["E"] > k1["Cminus"]) "as expected" else "AGAINST the expectation"))
p_null <- null_props(lengths(R$L_by_type))
rm_null <- rmse(matrix(p_null[types], n_spots, length(types), byrow = TRUE), truth)
score <- function(r) rmse(fill_failed(complete_prop(r$prop, rownames(truth), types), p_null), truth)
xP_real <- data$P; rownames(xP_real) <- rownames(R$mu_E)
k2 <- c(P_ref = score(run_nnls(data$P, RP$profiles)), real_ref = score(run_nnls(xP_real, R$profiles)))
cat(sprintf("K2 NNLS on P spots: RMSE with P reference %.4f, with real reference %.4f, null baseline %.4f -> expected P ref < null <= real ref: %s\n",
            k2["P_ref"], k2["real_ref"], rm_null, if (k2["P_ref"] < rm_null && k2["real_ref"] >= rm_null) "as expected" else "AGAINST the expectation"))
cat(sprintf("Time for arm P in one dataset (μ^P, φ^P, reference of %d cells, 500 spots): %.1f s\n", ncol(RP$ref), t_p))
if (any(res == "FAIL")) stop(sum(res == "FAIL"), " check(s) failed")
cat("All", length(res), "checks PASS\n")
