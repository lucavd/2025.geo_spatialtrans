# 03_reference.R — inputs of the generator and of the methods, per tissue (ADEMP §2.3, §2.6, §4)
# One pass over the count matrix of the raw h5ad (read in windows, R/stream_counts.R; the file is
# not changed, PROTOCOLLO §3.7). For each tissue it builds, on the kept genes:
#   - gene filter (ADEMP §2.6.2): count > 1 in >= 10 training cells, decided on training only;
#   - μ (DECISIONS 2026-10-29, noise model): per type, sum of the counts of gene g / sum of all
#     counts of those cells, on all genes; the kept genes are then taken, so a column sums to the
#     share of counts in kept genes (just under 1) and mean count = μ × L matches the real data;
#     μ_E from training (types >= 50 cells), μ_E′ from the E′ donor;
#   - L: library size (all genes) of every training cell, by type;
#   - φ_g: glmGamPoi on a sample of training cells, design ~ cell_type + donor, L as size factor
#     (DECISIONS 2026-10-29 noise model, donor covariate added by Luca 2026-10-30);
#   - method reference: up to ref_cap training cells per type (fixed seed), plus the mean counts
#     per cell of each type over all training cells (NNLS profiles);
#   - held-out cells (anchor) with their types.
# mode = "pilot" (ADEMP §10): held-out and E′ are drawn from the training donors with
# assign_split(qc = FALSE, seed 20261030); the main held-out and E′ donors are not read at all.
# mode = "main": the split of 02_split.R (held-out = first donor of the held-out group, E′ = second;
# the other group donors are used only in Aim 1).
#
# Pre-registered checks (written 2026-10-30 before the first run; the script stops if one fails):
#   1. pilot: no main held-out/E′ cell is used; pilot held-out, pilot E′ and pilot training are
#      disjoint by donor;
#   2. every type in the spots has >= 50 cells in training, held-out and E′ (ADEMP §2.6.1);
#   3. μ: no NA, columns in (0.9, 1];
#   4. φ: finite, >= 0, one per kept gene;
#   5. reference: each type has min(ref_cap, n training cells) cells; held-out cells all from one donor.
# Reported, not a stop: kept genes, φ range, library sizes, sums of n random training L
# (can a spot of 5 / 15 cells reach the depth targets 10,000 / 20,000?), reuse per held-out cell.
#
# Settings (Pierangela, 2026-10-30): method reference capped at 1,000 cells per type.
# φ: 20,000 training cells per tissue, stratified by type (DECISIONS 2026-10-29, refined 2026-10-30).
# Optional 3rd argument "all" fits φ on all training cells and saves <tissue>_<mode>_phiall.rds,
# used only by the pre-registered colon check (tests/test_phi_sample.R).
#
# Usage (working directory = thesis/):  Rscript scripts/03_reference.R <tissue> [pilot|main] [all]

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
mode <- if (length(args) >= 2) args[2] else "pilot"
stopifnot(tissue %in% c("brain", "colon", "lymphoid"), mode %in% c("pilot", "main"))

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
out_dir <- file.path(data_dir, "derived", "reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source("R/read_h5ad.R"); source("R/stream_counts.R"); source("R/assign_split.R")

min_cells <- 50          # ADEMP §2.6.1
gene_min_cells <- 10     # ADEMP §2.6.2: count > 1 in >= 10 training cells
ref_cap <- 1000          # cells per type in the method reference (Pierangela, 2026-10-30)
phi_all <- length(args) >= 3 && args[3] == "all"
phi_cells <- if (phi_all) Inf else 20000   # training cells for φ (DECISIONS 2026-10-30)
pilot_seed <- 20261030   # ADEMP §10
ref_seed <- 20261031     # reference and φ samples
counts_group <- c(brain = "/X", colon = "/raw/X", lymphoid = "/raw/X")[[tissue]]

t_start <- proc.time()[["elapsed"]]
cells <- read.csv(file.path(data_dir, "derived", "split", paste0(tissue, "_cells.csv.gz")),
                  stringsAsFactors = FALSE)

# Roles: train / heldout / eprime
if (mode == "pilot") {
  main_out <- unique(cells$donor_id[cells$split %in% c("heldout", "eprime")])
  cells <- cells[cells$split == "training", ]
  pd <- assign_split(cells, seed = pilot_seed, min_cells = min_cells, qc = FALSE)
  cells$role <- c(training = "train", heldout = "heldout", eprime = "eprime")[pd$split[match(cells$donor_id, pd$donor_id)]]
  stopifnot("Check 1: a main held-out/E′ donor in the pilot" = !any(cells$donor_id %in% main_out))
} else {
  cells$role <- c(training = "train", heldout = "heldout", eprime = "eprime")[cells$split]
}
stopifnot("Check 1: a donor in two roles" =
            all(tapply(cells$role, cells$donor_id, function(r) length(unique(r))) == 1))
donor_of <- function(role) unique(cells$donor_id[cells$role == role])
cat(sprintf("[%s, %s] training %d donors; held-out %s; E′ %s\n", tissue, mode,
            length(donor_of("train")), donor_of("heldout"), donor_of("eprime")))

# Cell types (>= 50 cells within each role)
n_by <- function(role) table(cells$cell_type[cells$role == role])
types_ref <- sort(names(which(n_by("train") >= min_cells)))
types_held <- sort(names(which(n_by("heldout") >= min_cells)))
types_ep <- sort(names(which(n_by("eprime") >= min_cells)))
types_spot <- Reduce(intersect, list(types_ref, types_held, types_ep))
cat(sprintf("  types: reference %d, held-out %d, E′ %d, spots %d (%s)\n", length(types_ref),
            length(types_held), length(types_ep), length(types_spot), paste(types_spot, collapse = ", ")))
cells <- cells[cells$cell_type %in% union(types_ref, union(types_held, types_ep)), ]

# Samples drawn before reading (fixed seed): reference cells and φ cells, from training
set.seed(ref_seed)
tr <- which(cells$role == "train" & cells$cell_type %in% types_ref)
ref_rows <- unlist(lapply(split(cells$row[tr], cells$cell_type[tr]), function(r)
  if (length(r) <= ref_cap) r else sample(r, ref_cap)))
phi_rows <- if (length(tr) <= phi_cells) cells$row[tr] else {
  # stratified by type, proportional to its training share
  unlist(lapply(split(cells$row[tr], cells$cell_type[tr]), function(r)
    sample(r, max(1, round(phi_cells * length(r) / length(tr))))))
}

# One pass over the counts
file <- list.files(file.path(data_dir, "raw"), pattern = paste0("^", tissue, "_.*\\.h5ad$"), full.names = TRUE)
genes <- read_var_genes(file, counts_group)
G <- nrow(genes)
key <- paste(cells$role, cells$cell_type)
groups <- sort(unique(key))
row_info <- cells[match(seq_len(max(cells$row)), cells$row), c("role", "cell_type")]
sums <- matrix(0, G, length(groups), dimnames = list(genes$gene_id, groups))
detect <- numeric(G)
L <- numeric(max(cells$row))   # library size by h5ad row
keep_rows <- unique(c(ref_rows, phi_rows, cells$row[cells$role == "heldout"]))
kept_mats <- list()
cat("  reading counts...\n")
stream_counts(file, counts_group, cells$row, G, function(m, r) {
  g <- match(paste(row_info$role[r], row_info$cell_type[r]), groups)
  ind <- Matrix::sparseMatrix(i = seq_along(r), j = g, x = 1, dims = c(length(r), length(groups)))
  sums <<- sums + as.matrix(Matrix::crossprod(m, ind))
  is_tr <- row_info$role[r] == "train"
  if (any(is_tr)) detect <<- detect + Matrix::colSums(m[is_tr, , drop = FALSE] > 1)
  L[r] <<- Matrix::rowSums(m)
  k <- r %in% keep_rows
  if (any(k)) { mk <- m[k, , drop = FALSE]; rownames(mk) <- r[k]; kept_mats[[length(kept_mats) + 1]] <<- mk }
})
kept <- do.call(rbind, kept_mats); rm(kept_mats); invisible(gc())
cat(sprintf("  read in %.1f min\n", (proc.time()[["elapsed"]] - t_start) / 60))

# Gene filter (training only)
keep_g <- detect >= gene_min_cells
cat(sprintf("  genes kept: %d of %d\n", sum(keep_g), G))
gk <- genes$gene_id[keep_g]

# μ and NNLS profiles
col_of <- function(role, types) paste(role, types)
tot <- function(role, types) vapply(types, function(k) sum(L[cells$row[cells$role == role & cells$cell_type == k]]), numeric(1))
mu_of <- function(role, types) {
  s <- sums[gk, col_of(role, types), drop = FALSE]
  out <- sweep(s, 2, tot(role, types), "/"); colnames(out) <- types; out
}
mu_E <- mu_of("train", types_ref)
mu_Ep <- mu_of("eprime", types_ep)
mu_held <- mu_of("heldout", types_held)   # Aim 1 target (not used by the pilot)
profiles <- sweep(sums[gk, col_of("train", types_ref), drop = FALSE], 2, as.numeric(n_by("train")[types_ref]), "/")
colnames(profiles) <- types_ref
for (nm in c("mu_E", "mu_Ep")) {
  cs <- colSums(get(nm))
  if (anyNA(get(nm)) || any(cs <= 0.9 | cs > 1 + 1e-9)) stop("Check 3: column sums of ", nm, " out of (0.9, 1]")
}
L_by_type <- lapply(setNames(types_ref, types_ref), function(k)
  unname(L[cells$row[cells$role == "train" & cells$cell_type == k]]))

# Method reference and held-out cells (genes x cells, kept genes)
take <- function(rows) { x <- Matrix::t(kept[as.character(rows), keep_g, drop = FALSE]); rownames(x) <- gk; x }
ref <- take(ref_rows)
ref_meta <- cells[match(ref_rows, cells$row), c("row", "donor_id", "cell_type")]
colnames(ref) <- paste0("c", ref_rows)
n_ref <- table(ref_meta$cell_type)
if (!all(n_ref[types_ref] == pmin(ref_cap, n_by("train")[types_ref]))) stop("Check 5: reference size per type")
held_rows <- cells$row[cells$role == "heldout"]
held <- take(held_rows); colnames(held) <- paste0("c", held_rows)
held_meta <- cells[match(held_rows, cells$row), c("row", "donor_id", "cell_type")]
if (length(unique(held_meta$donor_id)) != 1) stop("Check 5: held-out cells from more than one donor")

# φ_g with glmGamPoi (~ cell_type + donor, L as size factor)
t_phi <- proc.time()[["elapsed"]]
Y <- as.matrix(take(phi_rows))
cd <- data.frame(type = factor(cells$cell_type[match(phi_rows, cells$row)]),
                 donor = factor(cells$donor_id[match(phi_rows, cells$row)]))
fit <- glmGamPoi::glm_gp(Y, design = ~ type + donor, col_data = cd, size_factors = L[phi_rows])
phi <- setNames(fit$overdispersions, gk)
rm(Y, fit); invisible(gc())
if (length(phi) != length(gk) || any(!is.finite(phi)) || any(phi < 0)) stop("Check 4: φ")
cat(sprintf("  φ on %d cells in %.1f min; median %.2f, IQR %.2f–%.2f\n", length(phi_rows),
            (proc.time()[["elapsed"]] - t_phi) / 60, median(phi), quantile(phi, .25), quantile(phi, .75)))

# Reports for the factor levels (ADEMP §2.2) and the anchor reuse (§2.2 note)
set.seed(ref_seed + 1)
Ls <- unlist(L_by_type[types_spot])
for (n in c(5, 15)) {
  s <- replicate(10000, sum(sample(Ls, n, replace = TRUE)))
  cat(sprintf("  sum of %2d training library sizes (spot types): median %s; < 10,000: %.0f%%; < 20,000: %.0f%%\n",
              n, format(round(median(s)), big.mark = ","), 100 * mean(s < 1e4), 100 * mean(s < 2e4)))
}
held_spot <- sum(held_meta$cell_type %in% types_spot)
for (n in c(5, 15)) cat(sprintf("  anchor, 500 spots x %d cells: each held-out cell used %.1f times on average per replicate (%d cells)\n",
                                n, 500 * n / held_spot, held_spot))
cat(sprintf("  library size per cell (training): median %s, IQR %s–%s\n",
            format(round(median(L[cells$row[cells$role == "train"]])), big.mark = ","),
            format(round(quantile(L[cells$row[cells$role == "train"]], .25)), big.mark = ","),
            format(round(quantile(L[cells$row[cells$role == "train"]], .75)), big.mark = ",")))

res <- list(tissue = tissue, mode = mode, genes = genes[keep_g, ],
            types_ref = types_ref, types_held = types_held, types_ep = types_ep, types_spot = types_spot,
            mu_E = mu_E, mu_Ep = mu_Ep, mu_held = mu_held, profiles = profiles, L_by_type = L_by_type,
            phi = phi, ref = ref, ref_meta = ref_meta, held = held, held_meta = held_meta,
            donors = list(train = donor_of("train"), heldout = donor_of("heldout"), eprime = donor_of("eprime")),
            settings = list(min_cells = min_cells, gene_min_cells = gene_min_cells, ref_cap = ref_cap,
                            phi_cells = length(phi_rows), pilot_seed = pilot_seed, ref_seed = ref_seed),
            created = Sys.time())
saveRDS(res, file.path(out_dir, sprintf("%s_%s%s.rds", tissue, mode, if (phi_all) "_phiall" else "")))
cat(sprintf("  all checks passed; saved %s_%s.rds; total %.1f min\n", tissue, mode,
            (proc.time()[["elapsed"]] - t_start) / 60))
