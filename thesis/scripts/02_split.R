# 02_split.R — donor-level split per tissue (PROTOCOLLO §5.2, §5.6.3, §6 week 2; ADEMP §2.5, §2.6)
# Reads only the cell metadata (obs) and the row pointers of the count matrix from the raw
# h5ad files of 01_download.R; the files are not changed (§3.7). For each tissue it
#   1. keeps the biological subset (DECISIONS "Dataset", 2026-10-26, refined 2026-10-27):
#        brain    = PsychAD, disease "normal", age >= 18 years (unknown age excluded);
#        colon    = Pelka 2021, 10x 3' v2, MMRp tumours;
#        lymphoid = Wells 2025, thoracic lymph node, 10x 5' v2;
#   2. keeps cells with >= 100 detected genes and a non-missing author label (ADEMP §2.6.3);
#   3. assigns donors with the split rule (DECISIONS 2026-10-29, R/assign_split.R):
#      donor QC = donors with fewer cells than median - 3 MAD (log scale) are excluded;
#      held-out and E′ = two donors drawn at random (seed below) among the others;
#      all other donors -> training;
#   4. flags, within each split, the cell types with < 50 cells (ADEMP §2.6.1: not used for μ).
# Cell type = the authors' label (Pierangela, 2026-10-27): colon ClusterMidway,
# lymphoid author_cell_type, brain subclass.
# The gene filter (ADEMP §2.6.2) needs the counts and is decided on training only: it is done
# when μ is built, not here.
#
# Pre-registered checks (written before the first run, 2026-10-27; the script stops if one fails):
#   1. eligible donors: brain 208, colon 20, lymphoid 21; cells before the gene filter:
#      colon 89,242, lymphoid 185,437 (obs counts read on 2026-10-27/28);
#   2. every eligible donor has exactly one split; 1 held-out, 1 E′,
#      training = eligible - excluded - 2;
#   3. no donor_id in two splits (also in tests/test_split.R, PROTOCOLLO §5.6.3);
#   4. the count matrix has one row per obs cell.
# Amendment 2026-10-29, written before the run with the new rule (DECISIONS 2026-10-29): the
# 2026-10-27 coverage rule is replaced by donor QC + random draw, seed 20261029. The first run
# (2026-10-27) used the old rule; its donors are not used.
#   5. donors excluded by QC: brain 7, colon 0, lymphoid 0 (counted on the 2026-10-27 outputs).
# Amendment 2026-10-30, written before the run (DECISIONS 2026-10-30, held-out group): 4 donors are
# drawn instead of 2 (same seed); the 3rd and 4th complete the Aim 1 held-out group.
#   6. the first two drawn are unchanged: held-out/E′ = brain Donor_207/Donor_680,
#      colon C104/C105, lymphoid 583B/591C; training = eligible - excluded - 4.
# Expected, not a stop: in brain no cell is removed by the 100-genes filter (PsychAD minimum is
# 998 genes per nucleus in obs n_genes).
#
# Outputs in TESI_DATA/derived/split/ (outside git):
#   <tissue>_donors.csv    donor_id, n_cells, n_types_ge50, split (incl. "excluded")
#   <tissue>_cells.csv.gz  row (1-based row in the h5ad), cell_id, donor_id, sample_id,
#                          cell_type, n_genes, split, type_kept
#   <tissue>_types.csv     split, cell_type, n_cells, type_kept
#
# Run from the thesis project (working directory = thesis/):
#   source("scripts/02_split.R")       or   Rscript scripts/02_split.R

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
raw_dir <- file.path(data_dir, "raw")
out_dir <- file.path(data_dir, "derived", "split")
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
stopifnot("Package rhdf5 missing: renv::restore()" = requireNamespace("rhdf5", quietly = TRUE))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source("R/read_h5ad.R")
source("R/assign_split.R")

min_genes <- 100   # ADEMP §2.6.3
min_cells <- 50    # ADEMP §2.6.1
split_seed <- 20261029   # DECISIONS 2026-10-29; fixed before the run, never changed
n_draw <- 4              # held-out group (DECISIONS 2026-10-30)

# Age in years from the CELLxGENE development_stage label; month stages count as < 1 year,
# "unknown" becomes NA (excluded).
age_years <- function(stage) {
  age <- suppressWarnings(as.numeric(sub("-year.*", "", stage)))
  age[grepl("^80 year-old and over", stage)] <- 80
  age[grepl("-month-old stage$", stage)] <- 0
  age
}

tissues <- list(
  brain = list(
    counts = "/X",          # PsychAD: raw counts are in X
    label = "subclass", sample = NULL,
    cols = c("disease", "development_stage"),
    keep = function(o) o$disease == "normal" & !is.na(age_years(o$development_stage)) &
                       age_years(o$development_stage) >= 18,
    expected_donors = 208, expected_cells = NA, expected_excluded = 7,
    expected_pair = c("Donor_207", "Donor_680")),
  colon = list(
    counts = "/raw/X",      # X is normalised; raw counts in raw/X
    label = "ClusterMidway", sample = "biosample_id",
    cols = c("assay", "MMRStatusTumor", "disease"),
    keep = function(o) o$assay == "10x 3' v2" & o$MMRStatusTumor == "MMRp" &
                       o$disease == "colon adenocarcinoma",
    expected_donors = 20, expected_cells = 89242, expected_excluded = 0,
    expected_pair = c("C104", "C105")),
  lymphoid = list(
    counts = "/raw/X",
    label = "author_cell_type", sample = "library_id",
    cols = c("assay", "tissue"),
    keep = function(o) o$assay == "10x 5' v2" & o$tissue == "thoracic lymph node",
    expected_donors = 21, expected_cells = 185437, expected_excluded = 0,
    expected_pair = c("583B", "591C"))
)

summary_rows <- list()
for (tissue in names(tissues)) {
  cfg <- tissues[[tissue]]
  file <- list.files(raw_dir, pattern = paste0("^", tissue, "_.*\\.h5ad$"), full.names = TRUE)
  stopifnot("Exactly one raw h5ad per tissue" = length(file) == 1)
  cat(sprintf("\n[%s] %s\n", tissue, basename(file)))

  obs <- read_obs_columns(file, unique(c("donor_id", cfg$label, cfg$sample, cfg$cols)))
  obs$row <- seq_len(nrow(obs))
  n_genes <- n_genes_per_cell(file, cfg$counts)
  stopifnot("Check 4: count matrix rows != obs cells" = length(n_genes) == nrow(obs))
  obs$n_genes <- n_genes

  # 1. biological subset
  sel <- obs[cfg$keep(obs), ]
  n_donors <- length(unique(sel$donor_id))
  cat(sprintf("  subset: %d cells, %d donors\n", nrow(sel), n_donors))
  if (n_donors != cfg$expected_donors)
    stop(sprintf("Check 1: %d eligible donors, expected %d", n_donors, cfg$expected_donors))
  if (!is.na(cfg$expected_cells) && nrow(sel) != cfg$expected_cells)
    stop(sprintf("Check 1: %d cells, expected %d", nrow(sel), cfg$expected_cells))

  # 2. cell filters
  n_low <- sum(sel$n_genes < min_genes)
  n_nolab <- sum(is.na(sel[[cfg$label]]) | sel[[cfg$label]] == "")
  sel <- sel[sel$n_genes >= min_genes & !is.na(sel[[cfg$label]]) & sel[[cfg$label]] != "", ]
  cat(sprintf("  removed: %d cells with < %d genes, %d without label\n", n_low, min_genes, n_nolab))

  cells <- data.frame(
    row = sel$row, cell_id = sel$cell_id, donor_id = sel$donor_id,
    sample_id = if (is.null(cfg$sample)) sel$donor_id else sel[[cfg$sample]],  # PsychAD obs: no sample column, one DLPFC sample per donor_id
    cell_type = sel[[cfg$label]], n_genes = sel$n_genes, stringsAsFactors = FALSE)

  # 3. donor split
  donors <- assign_split(cells, seed = split_seed, min_cells = min_cells, n_draw = n_draw)
  n_excl <- sum(donors$split == "excluded")
  cat(sprintf("  donor QC: lower bound %.0f cells, %d donors excluded
",
              attr(donors, "qc_lower_bound"), n_excl))
  if (n_excl != cfg$expected_excluded)
    stop(sprintf("Check 5: %d donors excluded by QC, expected %d", n_excl, cfg$expected_excluded))
  stopifnot(
    "Check 2: donors lost or duplicated" = setequal(donors$donor_id, unique(cells$donor_id)),
    "Check 2: donor listed twice" = !anyDuplicated(donors$donor_id),
    "Check 2: not 1 held-out + 1 E′" = sum(donors$split == "heldout") == 1 && sum(donors$split == "eprime") == 1,
    "Check 6: training != eligible - excluded - 4" =
      sum(donors$split == "training") == n_donors - n_excl - n_draw,
    "Check 6: first two drawn donors changed" =
      identical(c(donors$donor_id[donors$split == "heldout"], donors$donor_id[donors$split == "eprime"]), cfg$expected_pair))
  cells$split <- donors$split[match(cells$donor_id, donors$donor_id)]
  cells <- cells[cells$split != "excluded", ]   # excluded donors stay only in <tissue>_donors.csv
  stopifnot("Check 3: a donor in two splits" =
              all(tapply(cells$split, cells$donor_id, function(s) length(unique(s))) == 1))

  # 4. cell types with < min_cells cells in a split
  types <- aggregate(list(n_cells = cells$row), cells[, c("split", "cell_type")], length)
  types$type_kept <- types$n_cells >= min_cells
  types <- types[order(types$split, -types$n_cells), ]
  cells$type_kept <- types$type_kept[match(paste(cells$split, cells$cell_type),
                                           paste(types$split, types$cell_type))]

  write.csv(donors, file.path(out_dir, paste0(tissue, "_donors.csv")), row.names = FALSE)
  write.csv(types, file.path(out_dir, paste0(tissue, "_types.csv")), row.names = FALSE)
  gz <- gzfile(file.path(out_dir, paste0(tissue, "_cells.csv.gz")), "w")
  write.csv(cells, gz, row.names = FALSE)
  close(gz)

  for (s in c("training", "heldout", "eprime", paste0("group", 3:n_draw))) {
    tt <- types[types$split == s, ]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      tissue = tissue, split = s,
      donors = if (s == "training") sum(donors$split == s) else donors$donor_id[donors$split == s],
      cells = sum(cells$split == s), types_kept = sum(tt$type_kept), types_dropped = sum(!tt$type_kept))
  }
}

cat("\nSplit summary (cell types kept = >= 50 cells in that split):\n")
print(do.call(rbind, summary_rows), row.names = FALSE)
cat("\nAll pre-registered checks passed. Outputs in", normalizePath(out_dir), "\n")
