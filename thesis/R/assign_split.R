# assign_split.R — donor-level split rule (DECISIONS 2026-10-29, ADEMP §2.5).
# 1. Donor QC: a donor is excluded when its number of cells is a low outlier, i.e. below
#    median - 3 MAD on the log scale (the scater::isOutlier convention for cells, OSCA,
#    Amezquita et al. 2020, doi:10.1038/s41592-019-0654-x, applied here to donors).
#    Excluded donors go to no split.
# 2. Among the remaining donors, two are drawn at random with a fixed seed: the first is the
#    held-out donor, the second E′. All other donors go to training.
# Donors are sorted by donor_id before the draw, so the result depends only on the seed and on
# the set of donors, not on the order of the cells.

# cells: data.frame with one row per cell and columns donor_id, cell_type.
# Returns one row per donor: donor_id, n_cells, n_types_ge<min_cells>, split
# (training, heldout, eprime or excluded).
# qc = FALSE skips step 1: used for the pilot split, drawn from training donors that already
# passed QC (ADEMP §10).
assign_split <- function(cells, seed, min_cells = 50, n_mads = 3, qc = TRUE) {
  stopifnot(all(c("donor_id", "cell_type") %in% names(cells)), !anyNA(cells$donor_id),
            "seed must be one number" = length(seed) == 1 && is.numeric(seed))
  per_donor <- split(cells$cell_type, cells$donor_id)
  donors <- data.frame(
    donor_id = names(per_donor),
    n_cells = vapply(per_donor, length, integer(1)),
    n_types = vapply(per_donor, function(x) sum(table(x) >= min_cells), integer(1)),
    stringsAsFactors = FALSE
  )
  donors <- donors[order(donors$donor_id), ]

  log_cells <- log(donors$n_cells)
  lower <- if (qc) exp(median(log_cells) - n_mads * mad(log_cells)) else 0
  donors$split <- ifelse(donors$n_cells < lower, "excluded", "training")
  kept <- which(donors$split == "training")
  stopifnot("At least 3 donors needed after QC (training, held-out, E′)" = length(kept) >= 3)

  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")   # same draw in every R version >= 3.6
  set.seed(seed)
  drawn <- kept[sample.int(length(kept), 2)]
  donors$split[drawn] <- c("heldout", "eprime")

  names(donors)[names(donors) == "n_types"] <- paste0("n_types_ge", min_cells)
  attr(donors, "qc_lower_bound") <- lower
  rownames(donors) <- NULL
  donors
}
