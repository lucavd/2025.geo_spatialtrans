# 07_aim1.R — Aim 1, fidelity of the cell-type profiles μ of each source to real held-out profiles
# (PROTOCOLLO §1, §5.5; ADEMP §5.1, §9; DECISIONS 2026-10-26, 2026-10-30, 2026-11-04, 2026-11-11).
#
# Design (fixed by the ADEMP; implementation choices written 2026-11-11 before the first run):
#   targets: each of the 4 donors of the held-out group (DECISIONS 2026-10-30); per target, the types with
#     >= 50 cells in that donor and in the training reference (ADEMP §2.6.1; others "not evaluated");
#   arms: E (training μ), E′ (each other donor of the group as source: 12 ordered pairs = 6 pairs, the
#     correlations and Jaccard are symmetric, the silhouette is not), C− (μ_E with gene labels permuted
#     within type, a new permutation per draw), P (de novo world, distribution properties only);
#     O is not an Aim 1 arm (its profile is the held-out profile by construction, ADEMP §2.1);
#     G only after its DECISIONS row;
#   per arm x target x type: Pearson (main) and Spearman of log(1 + 10^4 μ) over the kept genes;
#     Jaccard of the top-k markers (k = 10, 20, 50, 100, 200; main k = 100), markers = highest log2
#     fold-change vs the mean of the other types of the same comparison (the type set is the same in
#     source and target); silhouette width per type of cells simulated from the source with the common
#     noise model (as many per type as the target donor has, library sizes drawn from the training cells
#     of the type), against the silhouette of the real cells of the target donor on the same types;
#   P: separability (silhouette) and distributions of the gene means (log-profile) and of the log2
#     fold-changes of the top-100 markers per type: quantiles 10/50/90/99/99.9% and two-sample KS
#     distance, against the target donor and against E;
#   noise ceiling per target x type: split halves, Spearman-Brown, sqrt (Pearson and Spearman);
#   number of draws (DECISIONS 2026-11-11): B of the ceiling, n of C−, n of P = max(10,
#     ceiling((SD_upper / 0.005)^2)), SD_upper = upper 95% limit of the SD over a pilot of 10 draws,
#     the largest over donors, types and measures; computed in this tissue.
# Implementation choices (Pierangela, 2026-11-11, before any result):
#   - HVGs: Seurat VST reimplemented on streamed sums (R/aim1.R; check against Seurat in
#     tests/test_aim1.R), 2,000 genes among the kept genes (§2.6 list), on the raw counts of the
#     training cells of the reference types; the same HVGs for the held-out, E, E′ and C−;
#   - P has no real genes, so its HVGs are chosen with the same VST on its simulated method reference
#     (its "training", ADEMP §2.1), in each draw;
#   - measures that enter the rule of the number of draws: C− Pearson, Spearman, Jaccard k = 100,
#     silhouette; P KS distances (means and fold-changes, vs target and vs E) and silhouette; the
#     quantiles are descriptive (they are on the log-profile scale, not on the 0–1 scale of 0.005);
#   - a ceiling with negative Spearman-Brown reliability is set to 0 (not expected with >= 50 cells);
#   - streams: L'Ecuyer-CMRG from the study seed 42 (ADEMP §8), one stream per split / draw, so a
#     result does not depend on how many draws come before it in another component.
#
# Modes. main: the held-out group of 02_split.R. Blocked until the ADEMP is approved (as 04_simulate.R),
#   because Aim 1 on the held-out donors is a result. pilot (code test, ADEMP §10): 4 donors drawn from
#   the training donors with assign_split(qc = FALSE, seed 20261030, n_draw = 4) (the first two are the
#   pilot held-out and E′ of 03_reference.R), the main group is never read; the training is the main
#   training minus these 4 donors; φ, the kept genes and the reference of P come from <tissue>_pilot.rds
#   and the P calibration from the main training (both include pilot group donors: declared, test only).
#
# Pre-registered checks (the script stops if one fails):
#   A1 4 distinct group donors, none among the training cells read; the first two are the held-out and E′
#      donors of the reference file;
#   A2 μ of the training and of every group donor: no NA, columns in (0.9, 1];
#   A3 main only: μ_E, μ of donor 1 and of donor 2 recomputed here equal mu_E, mu_held, mu_Ep of
#      03_reference.R (max abs difference < 1e-12); always: same genes in the same order;
#   A4 a profile against itself gives Pearson = Spearman = Jaccard = 1;
#   A5 2,000 distinct HVGs, all among the kept genes;
#   A6 every ceiling in [0, 1];
#   A7 every number of draws < 3,000 (stream layout).
# Pre-registered expectations (reported, not a stop; ADEMP §1, §9):
#   H1: C− has the lowest Pearson among E, E′, C− for every target x type (C− Pearson near 0);
#   §9.2: E or E′ above the ceiling by more than 2 x MCSE of the ceiling is flagged for investigation.
#
# Usage (working directory = thesis/):  Rscript scripts/07_aim1.R <tissue> [pilot|main]
# Output: <TESI_DATA>/derived/aim1/<tissue>_<mode>_aim1.rds and CSVs (fidelity, silhouette, ceiling,
# p_dist, draws) next to it.

ADEMP_APPROVED <- NA_character_

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
mode <- if (length(args) >= 2) args[2] else "pilot"
run <- if (length(args) >= 3) args[3] else "run"        # plan = only the 10-draw pilots and the numbers of draws
workers <- if (length(args) >= 4) as.integer(args[4]) else 1L
stopifnot(tissue %in% c("colon", "lymphoid", "brain"), mode %in% c("pilot", "main"), run %in% c("plan", "run"),
          workers >= 1)
# Draws in parallel by fork (Linux server); on Windows they run one after the other. Every split / draw sets
# its own stream, so the results do not depend on the number of workers.
par_lapply <- function(X, FUN) {
  if (workers > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(X, FUN, mc.cores = workers, mc.preschedule = FALSE) else lapply(X, FUN)
}
if (mode == "main" && is.na(ADEMP_APPROVED)) stop("Main run blocked: the ADEMP is not approved yet")

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
source("renv/activate.R")
suppressPackageStartupMessages({ library(Matrix); library(scater) })
for (f in c("R/read_h5ad.R", "R/stream_counts.R", "R/assign_split.R", "R/simulate_spots.R",
            "R/arm_p.R", "R/aim1.R")) source(f)
out_dir <- file.path(data_dir, "derived", "aim1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]
elapsed <- function() sprintf("%.1f min", (proc.time()[["elapsed"]] - t0) / 60)

min_cells <- 50; pilot_seed <- 20261030; study_seed <- 42
ks <- c(10, 20, 50, 100, 200); k_main <- 100; n_hvg <- 2000; mcse_target <- 0.005; n_pilot <- 10
counts_group <- c(brain = "/X", colon = "/raw/X", lymphoid = "/raw/X")[[tissue]]

R0 <- readRDS(file.path(data_dir, "derived", "reference", sprintf("%s_%s.rds", tissue, mode)))
calib <- readRDS(file.path(data_dir, "derived", "p_calibration", sprintf("%s_p_calibration.rds", tissue)))
genes_k <- rownames(R0$mu_E); types_ref <- R0$types_ref
stopifnot(all(types_ref %in% names(calib$markers)))

# ---- Donors and cells
cells <- read.csv(file.path(data_dir, "derived", "split", paste0(tissue, "_cells.csv.gz")), stringsAsFactors = FALSE)
if (mode == "main") {
  grp <- c("heldout", "eprime", "group3", "group4")
  donors <- vapply(grp, function(s) unique(cells$donor_id[cells$split == s]), "")
  train_rows <- cells$row[cells$split == "training"]
} else {
  tr <- cells[cells$split == "training", ]
  pd <- assign_split(tr, seed = pilot_seed, min_cells = min_cells, qc = FALSE, n_draw = 4)
  donors <- vapply(c("heldout", "eprime", "group3", "group4"), function(s) pd$donor_id[pd$split == s], "")
  train_rows <- tr$row[!tr$donor_id %in% donors]
}
n_row <- max(cells$row)
cells <- cells[match(seq_len(n_row), cells$row), ]   # index by h5ad row (NA where no cell)
train_rows <- train_rows[cells$cell_type[train_rows] %in% types_ref]
group_rows <- which(cells$donor_id %in% donors)
if (length(unique(donors)) != 4 || any(cells$donor_id[train_rows] %in% donors) ||
    !identical(unname(donors[1:2]), c(R0$donors$heldout, R0$donors$eprime)))
  stop("A1: group donors")
cat(sprintf("[%s, %s] group %s; training %d cells, %d donors\n", tissue, mode, paste(donors, collapse = ", "),
            length(train_rows), length(unique(cells$donor_id[train_rows]))))

# ---- Pass 1: training sums for μ_E and VST; group donor cells kept (kept genes)
file <- list.files(file.path(data_dir, "raw"), pattern = paste0("^", tissue, "_.*\\.h5ad$"), full.names = TRUE)
genes_all <- read_var_genes(file, counts_group)
gi <- match(genes_k, genes_all$gene_id)
stopifnot("A3: kept genes not in the h5ad" = !anyNA(gi))
G <- length(genes_k)
is_train <- logical(n_row); is_train[train_rows] <- TRUE
s1 <- numeric(G); s2 <- numeric(G)
sum_type <- matrix(0, G, length(types_ref), dimnames = list(genes_k, types_ref))
L <- numeric(n_row)
grp_mats <- list()
cat("  pass 1...\n")
stream_counts(file, counts_group, c(train_rows, group_rows), nrow(genes_all), function(m, r) {
  L[r] <<- Matrix::rowSums(m)
  mk <- m[, gi, drop = FALSE]
  t <- is_train[r]
  if (any(t)) {
    mt <- mk[t, , drop = FALSE]
    s1 <<- s1 + Matrix::colSums(mt); s2 <<- s2 + Matrix::colSums(mt^2)
    j <- match(cells$cell_type[r[t]], types_ref)
    ind <- Matrix::sparseMatrix(i = seq_along(j), j = j, x = 1, dims = c(length(j), length(types_ref)))
    sum_type <<- sum_type + as.matrix(Matrix::crossprod(mt, ind))
  }
  if (any(!t)) { g <- mk[!t, , drop = FALSE]; rownames(g) <- r[!t]; grp_mats[[length(grp_mats) + 1]] <<- g }
}, verbose = FALSE)
cat("  pass 1 done,", elapsed(), "\n")

n_tr <- length(train_rows)
L_tr <- split(L[train_rows], cells$cell_type[train_rows])[types_ref]
mu_E <- sweep(sum_type, 2, vapply(L_tr, sum, 0), "/")
L_by_type <- lapply(L_tr, unname)
grp <- do.call(rbind, grp_mats); rm(grp_mats)
grp_row <- as.integer(rownames(grp))
X <- lapply(setNames(donors, donors), function(d) {
  i <- which(cells$donor_id[grp_row] == d)
  x <- Matrix::t(grp[i, , drop = FALSE]); rownames(x) <- genes_k
  ty <- cells$cell_type[grp_row[i]]
  n <- table(ty)
  types <- sort(intersect(names(n)[n >= min_cells], types_ref))
  list(x = x, type = ty, L = L[grp_row[i]], types = types,
       mu = pooled_mu(x, ty, L[grp_row[i]], sort(names(n)[n >= min_cells])))
})
rm(grp); invisible(gc())
for (m in c(list(mu_E), lapply(X, `[[`, "mu"))) {
  cs <- colSums(m)
  if (anyNA(m) || any(cs <= 0.9 | cs > 1 + 1e-9)) stop("A2: μ column sums out of (0.9, 1]")
}
if (mode == "main") {
  dif <- c(max(abs(mu_E - R0$mu_E[, types_ref])),
           max(abs(X[[1]]$mu[, R0$types_held] - R0$mu_held)), max(abs(X[[2]]$mu[, R0$types_ep] - R0$mu_Ep)))
  if (any(dif > 1e-12)) stop("A3: μ differs from 03_reference.R: ", paste(signif(dif, 3), collapse = ", "))
}
for (d in donors) cat(sprintf("  %s: %d cells, %d evaluated types (%s)\n", d, ncol(X[[d]]$x),
                              length(X[[d]]$types), paste(X[[d]]$types, collapse = ", ")))

# ---- Pass 2: VST standardized variance on the training cells
mean_g <- s1 / n_tr
var_g <- (s2 - n_tr * mean_g^2) / (n_tr - 1)
sd_g <- sqrt(vst_expected(mean_g, var_g))
sz <- numeric(G); nnz <- numeric(G)
cat("  pass 2 (VST)...\n")
stream_counts(file, counts_group, train_rows, nrow(genes_all), function(m, r) {
  p <- vst_pass2(m[, gi, drop = FALSE], mean_g, sd_g, sqrt(n_tr))
  sz <<- sz + p$sumz2; nnz <<- nnz + p$nnz
}, verbose = FALSE)
V <- vst_finish(mean_g, sd_g, n_tr, sz, nnz, n_hvg, genes_k)
hvg <- V$hvg
if (length(unique(hvg)) != n_hvg || !all(hvg %in% genes_k)) stop("A5: HVGs")
cat("  pass 2 done,", elapsed(), "\n")

# ---- Streams
RNGkind("L'Ecuyer-CMRG")
set.seed(study_seed)
base <- .Random.seed
for (k in seq_len(50000 + match(tissue, c("colon", "lymphoid", "brain")) * 10000)) base <- parallel::nextRNGStream(base)
stream_at <- function(k) { s <- base; for (j in seq_len(k)) s <- parallel::nextRNGStream(s); s }
use_stream <- function(k) assign(".Random.seed", stream_at(k), envir = .GlobalEnv)
off <- c(ceiling = 0, fixed = 3000, Cminus = 4000, P = 7000)

# ---- Helpers
n_target <- function(d, S) { n <- table(X[[d]]$type)[S]; setNames(as.integer(n), S) }
held_sil_cache <- list()
held_sil <- function(d, S) {
  key <- paste(d, paste(S, collapse = "|"))
  if (is.null(held_sil_cache[[key]])) {
    i <- X[[d]]$type %in% S
    held_sil_cache[[key]] <<- silhouette_by_type(X[[d]]$x[, i], X[[d]]$type[i], hvg)
  }
  held_sil_cache[[key]]
}
score <- function(src, d, arm, source, draw, phi = R0$phi, hvg_src = hvg) {
  S <- intersect(X[[d]]$types, colnames(src))
  tg <- X[[d]]$mu[, S, drop = FALSE]; sr <- src[, S, drop = FALSE]
  f <- merge(fidelity_cor(sr, tg), jaccard_topk(sr, tg, ks), by = "cell_type")
  sim <- simulate_type_cells(sr, phi, n_target(d, S), L_by_type)
  sil <- silhouette_by_type(sim$x, sim$type, hvg_src)
  f$silhouette <- sil[f$cell_type]
  f$silhouette_target <- held_sil(d, S)[f$cell_type]
  data.frame(target = d, arm = arm, source = source, draw = draw, n_types = length(S), f)
}

# A4: a profile against itself
d1 <- donors[[1]]
self <- merge(fidelity_cor(X[[d1]]$mu, X[[d1]]$mu), jaccard_topk(X[[d1]]$mu, X[[d1]]$mu, ks), by = "cell_type")
if (any(abs(as.matrix(self[, -1]) - 1) > 1e-12)) stop("A4: self-comparison not 1")

# ---- Noise ceiling (pilot of 10 splits, then up to B)
ceil_draw <- function(b) {
  use_stream(off[["ceiling"]] + b)
  do.call(rbind, lapply(donors, function(d) data.frame(target = d, split = b,
    noise_ceiling_draw(X[[d]]$x, X[[d]]$type, X[[d]]$L, X[[d]]$types))))
}
# SD_upper per measure x target x type over the pilot draws, and the rows that set n (largest SD_upper)
sd_table <- function(df, cols) {
  out <- do.call(rbind, lapply(cols, function(v) {
    s <- tapply(df[[v]], paste(df$target, df$cell_type, sep = " / "), sd_upper)
    data.frame(measure = v, target_type = names(s), sd = as.numeric(tapply(df[[v]], paste(df$target, df$cell_type, sep = " / "), sd)),
               sd_upper = as.numeric(s), n_rule = ceiling(as.numeric(s)^2 / mcse_target^2))
  }))
  out[order(-out$sd_upper), ]
}
show_driver <- function(tab, what) {
  cat(sprintf("  %s: largest SD_upper per measure, and the 5 target x type that set n:\n", what))
  print(aggregate(cbind(sd_upper, n_rule) ~ measure, tab, max), row.names = FALSE, digits = 3)
  print(head(tab, 5), row.names = FALSE, digits = 3)
}
t_ceil <- proc.time()[["elapsed"]]
ceil <- do.call(rbind, par_lapply(seq_len(n_pilot), ceil_draw))
min_per <- list(ceiling = (proc.time()[["elapsed"]] - t_ceil) / 60 / n_pilot)
sd_c_tab <- sd_table(ceil, c("ceiling_pearson", "ceiling_spearman"))
sd_ceil <- sd_c_tab$sd_upper
B <- n_draws_rule(sd_ceil, mcse_target, n_pilot)
if (B >= 3000) stop("A7: B")
show_driver(sd_c_tab, sprintf("noise ceiling, B = %d", B))
if (B > n_pilot && run == "run") ceil <- rbind(ceil, do.call(rbind, par_lapply((n_pilot + 1):B, ceil_draw)))
if (any(ceil$ceiling_pearson < 0 | ceil$ceiling_pearson > 1 | ceil$ceiling_spearman < 0 |
        ceil$ceiling_spearman > 1)) stop("A6: ceiling out of [0, 1]")
ceil_sum <- do.call(rbind, lapply(split(ceil, paste(ceil$target, ceil$cell_type)), function(z) data.frame(
  target = z$target[1], cell_type = z$cell_type[1], n_cells = sum(X[[z$target[1]]]$type == z$cell_type[1]),
  ceiling_pearson = mean(z$ceiling_pearson), mcse_pearson = sd(z$ceiling_pearson) / sqrt(nrow(z)),
  ceiling_spearman = mean(z$ceiling_spearman), mcse_spearman = sd(z$ceiling_spearman) / sqrt(nrow(z)), B = nrow(z))))
cat(sprintf("  noise ceiling: B = %d (SD_upper max %.4f), %s\n", B, max(sd_ceil), elapsed()))

# Silhouettes of the real target cells on each target's own types, computed once before the parallel draws
# (a cache filled inside a forked worker would be lost)
for (d in donors) invisible(held_sil(d, X[[d]]$types))

# ---- Fixed arms: E and E′ (one simulation each for the silhouette); skipped by plan
use_stream(off[["fixed"]] + 1)
fid <- list()
for (d in if (run == "run") donors else character(0)) {
  fid[[length(fid) + 1]] <- score(mu_E, d, "E", "training", 0)
  for (s in setdiff(donors, d)) fid[[length(fid) + 1]] <- score(X[[s]]$mu, d, "Eprime", s, 0)
}
cat("  E and E′ done,", elapsed(), "\n")

# ---- C− draws
cminus_draw <- function(i) {
  use_stream(off[["Cminus"]] + i)
  mu_C <- permute_mu(mu_E)
  do.call(rbind, lapply(donors, function(d) score(mu_C, d, "Cminus", "training", i)))
}
# ---- P draws (world as in 04_simulate.R: μ^P, φ^P by rank, simulated reference)
mean_real <- pooled_mean(R0$mu_E, R0$L_by_type)
L_ref <- Matrix::colSums(R0$ref)
p_draw <- function(i) {
  use_stream(off[["P"]] + i)
  P <- simulate_mu_p(calib, G, types_ref)
  phi_p <- phi_by_rank(R0$phi, mean_real, P$base)
  refP <- simulate_reference_p(P$mu, phi_p, R0$ref_meta, L_ref)
  hvg_p <- vst_dense(refP$ref, n_hvg)$hvg
  rm(refP)
  lfc_top <- function(mu, S) { l <- lfc_matrix(mu[, S, drop = FALSE]); lapply(setNames(seq_along(S), S), function(j) sort(l[, j], decreasing = TRUE)[seq_len(k_main)]) }
  out_s <- list(); out_d <- list()
  for (d in donors) {
    S <- X[[d]]$types
    sim <- simulate_type_cells(P$mu[, S, drop = FALSE], phi_p, n_target(d, S), L_by_type)
    sil <- silhouette_by_type(sim$x, sim$type, hvg_p)
    out_s[[d]] <- data.frame(target = d, arm = "P", source = "P", draw = i, n_types = length(S), cell_type = S,
                             silhouette = sil[S], silhouette_target = held_sil(d, S)[S])
    lp_P <- log_profile(P$mu); lp_T <- log_profile(X[[d]]$mu); lp_E <- log_profile(mu_E)
    fP <- lfc_top(P$mu, S); fT <- lfc_top(X[[d]]$mu, S); fE <- lfc_top(mu_E, S)
    out_d[[d]] <- do.call(rbind, lapply(S, function(k) data.frame(
      target = d, draw = i, cell_type = k,
      ks_mean_target = ks_d(lp_P[, k], lp_T[, k]), ks_mean_E = ks_d(lp_P[, k], lp_E[, k]),
      ks_fc_target = ks_d(fP[[k]], fT[[k]]), ks_fc_E = ks_d(fP[[k]], fE[[k]]),
      t(setNames(dist_quantiles(lp_P[, k]), paste0("mean_P_", names(dist_quantiles(1))))),
      t(setNames(dist_quantiles(lp_T[, k]), paste0("mean_target_", names(dist_quantiles(1))))),
      t(setNames(dist_quantiles(lp_E[, k]), paste0("mean_E_", names(dist_quantiles(1))))),
      t(setNames(dist_quantiles(fP[[k]]), paste0("fc_P_", names(dist_quantiles(1))))),
      t(setNames(dist_quantiles(fT[[k]]), paste0("fc_target_", names(dist_quantiles(1))))),
      t(setNames(dist_quantiles(fE[[k]]), paste0("fc_E_", names(dist_quantiles(1))))))))
  }
  list(sil = do.call(rbind, out_s), dist = do.call(rbind, out_d))
}

# ---- Pilots of 10 draws (C−, P), saved before anything else, then the numbers of draws
stem <- file.path(out_dir, sprintf("%s_%s_aim1", tissue, mode))
t_c <- proc.time()[["elapsed"]]
cm <- do.call(rbind, par_lapply(seq_len(n_pilot), cminus_draw))
min_per$Cminus <- (proc.time()[["elapsed"]] - t_c) / 60 / n_pilot
t_p <- proc.time()[["elapsed"]]
pp <- par_lapply(seq_len(n_pilot), p_draw)
min_per$P <- (proc.time()[["elapsed"]] - t_p) / 60 / n_pilot
p_sil <- do.call(rbind, lapply(pp, `[[`, "sil")); p_dist <- do.call(rbind, lapply(pp, `[[`, "dist"))
saveRDS(list(ceiling = ceil, cminus = cm, p_silhouette = p_sil, p_dist = p_dist, donors = donors, created = Sys.time()),
        paste0(stem, "_pilot_draws.rds"))
sd_C_tab <- sd_table(cm, c("pearson", "spearman", paste0("jaccard_k", k_main), "silhouette"))
sd_P_tab <- rbind(sd_table(p_sil, "silhouette"), sd_table(p_dist, c("ks_mean_target", "ks_mean_E", "ks_fc_target", "ks_fc_E")))
sd_P_tab <- sd_P_tab[order(-sd_P_tab$sd_upper), ]
sd_C <- sd_C_tab$sd_upper; sd_P <- sd_P_tab$sd_upper
n_C <- n_draws_rule(sd_C, mcse_target, n_pilot); n_P <- n_draws_rule(sd_P, mcse_target, n_pilot)
if (n_C >= 3000) stop("A7: n of C−")
if (n_P >= 3000) stop("A7: n of P")
show_driver(sd_C_tab, sprintf("C−, n = %d", n_C))
show_driver(sd_P_tab, sprintf("P, n = %d", n_P))
plan <- data.frame(component = c("ceiling", "Cminus", "P"), n = c(B, n_C, n_P),
                   wall_min_per_draw = unlist(min_per[c("ceiling", "Cminus", "P")]),
                   workers = if (.Platform$OS.type == "unix") workers else 1L)
plan$core_h_remaining <- plan$wall_min_per_draw * plan$workers * pmax(plan$n - n_pilot, 0) / 60
cat(sprintf("\nPlan (DECISIONS 2026-11-11), %s; times measured on the pilot draws:\n", elapsed()))
print(plan, row.names = FALSE, digits = 3)
write.csv(plan, paste0(stem, "_plan.csv"), row.names = FALSE)
if (run == "plan") { cat("Plan only: pilot draws in", paste0(stem, "_pilot_draws.rds"), "\n"); quit(save = "no") }

if (n_C > n_pilot) cm <- rbind(cm, do.call(rbind, par_lapply((n_pilot + 1):n_C, cminus_draw)))
cat("  C− done,", elapsed(), "\n")
if (n_P > n_pilot) {
  pp <- par_lapply((n_pilot + 1):n_P, p_draw)
  p_sil <- rbind(p_sil, do.call(rbind, lapply(pp, `[[`, "sil")))
  p_dist <- rbind(p_dist, do.call(rbind, lapply(pp, `[[`, "dist")))
}
cat("  P done,", elapsed(), "\n")

# ---- Tables
fid <- rbind(do.call(rbind, fid), cm)
draws <- data.frame(component = c("ceiling", "Cminus", "P"), n = c(B, n_C, n_P),
                    sd_upper_max = c(max(sd_ceil), max(sd_C, na.rm = TRUE), max(sd_P, na.rm = TRUE)))
# Summary per arm x target x type (mean over draws, MCSE for C−)
agg <- do.call(rbind, lapply(split(fid, paste(fid$arm, fid$target, fid$cell_type, fid$source)), function(z) {
  data.frame(arm = z$arm[1], target = z$target[1], source = z$source[1], cell_type = z$cell_type[1], n_draws = nrow(z),
             pearson = mean(z$pearson), mcse_pearson = if (nrow(z) > 1) sd(z$pearson) / sqrt(nrow(z)) else NA,
             spearman = mean(z$spearman), jaccard_k100 = mean(z$jaccard_k100),
             silhouette = mean(z$silhouette), silhouette_target = z$silhouette_target[1])
}))
rownames(agg) <- NULL
res <- list(tissue = tissue, mode = mode, donors = donors, hvg = hvg, types = lapply(X, `[[`, "types"),
            fidelity = fid, summary = agg, ceiling = ceil, ceiling_summary = ceil_sum, p_silhouette = p_sil,
            p_dist = p_dist, draws = draws, settings = list(seed = study_seed, ks = ks, n_hvg = n_hvg,
            mcse_target = mcse_target, n_pilot = n_pilot, min_cells = min_cells), created = Sys.time())
stem <- file.path(out_dir, sprintf("%s_%s_aim1", tissue, mode))
saveRDS(res, paste0(stem, ".rds"))
write.csv(fid, paste0(stem, "_fidelity.csv"), row.names = FALSE)
write.csv(agg, paste0(stem, "_summary.csv"), row.names = FALSE)
write.csv(ceil_sum, paste0(stem, "_ceiling.csv"), row.names = FALSE)
write.csv(p_sil, paste0(stem, "_p_silhouette.csv"), row.names = FALSE)
write.csv(p_dist, paste0(stem, "_p_dist.csv"), row.names = FALSE)
write.csv(draws, paste0(stem, "_draws.csv"), row.names = FALSE)

# ---- Report
cat("\nNumber of draws (DECISIONS 2026-11-11):\n"); print(draws, row.names = FALSE)
w <- reshape(agg[agg$arm != "Eprime", c("arm", "target", "cell_type", "pearson")], idvar = c("target", "cell_type"),
             timevar = "arm", direction = "wide")
ep <- aggregate(pearson ~ target + cell_type, agg[agg$arm == "Eprime", ], mean)
names(ep)[3] <- "pearson.Eprime_mean"
w <- merge(merge(w, ep, all.x = TRUE), ceil_sum[, c("target", "cell_type", "n_cells", "ceiling_pearson", "mcse_pearson")])
cat("\nPearson on log-profiles per target x type (E′ = mean over the other group donors):\n")
print(format(w, digits = 3), row.names = FALSE)

# H1 and §9.2
cmp <- merge(agg[agg$arm == "Cminus", c("target", "cell_type", "pearson")],
             aggregate(pearson ~ target + cell_type, agg[agg$arm %in% c("E", "Eprime"), ], min), by = c("target", "cell_type"),
             suffixes = c("_C", "_min_real"))
cat(sprintf("\nH1: C− lowest Pearson in %d of %d target x type (C− Pearson range %.3f to %.3f)\n",
            sum(cmp$pearson_C < cmp$pearson_min_real), nrow(cmp), min(cmp$pearson_C), max(cmp$pearson_C)))
over <- merge(agg[agg$arm %in% c("E", "Eprime"), c("arm", "target", "source", "cell_type", "pearson")],
              ceil_sum[, c("target", "cell_type", "ceiling_pearson", "mcse_pearson")])
over <- over[over$pearson > over$ceiling_pearson + 2 * over$mcse_pearson, ]
cat(sprintf("§9.2: %d comparisons of E/E′ above the ceiling + 2 MCSE%s\n", nrow(over),
            if (nrow(over)) " (to investigate, listed below)" else ""))
if (nrow(over)) print(over, row.names = FALSE)

s <- aggregate(cbind(jaccard_k100, silhouette, silhouette_target) ~ arm, agg, mean)
ps <- aggregate(cbind(silhouette, silhouette_target) ~ arm, p_sil, mean); ps$jaccard_k100 <- NA
cat("\nMean over targets x types: Jaccard top-100, silhouette of the arm, silhouette of the real target cells:\n")
print(format(rbind(s, ps[, names(s)]), digits = 3), row.names = FALSE)
pd_s <- colMeans(p_dist[, c("ks_mean_target", "ks_mean_E", "ks_fc_target", "ks_fc_E")])
cat("\nP, mean KS distance:", paste(names(pd_s), sprintf("%.3f", pd_s), collapse = "; "), "\n")
cat(sprintf("\nWritten %s.rds and CSVs; total %s\n", stem, elapsed()))
