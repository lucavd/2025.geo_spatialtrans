# test_split.R — tests of the donor split (PROTOCOLLO §5.6.3: no donor in two splits).
# Part 1 tests assign_split() on toy data; part 2 checks the real outputs of 02_split.R.
# Run from the thesis project:  source("tests/test_split.R")

source("R/assign_split.R")

# Part 1: toy data with known answer.
# Six donors with 900-1200 cells and one tiny donor (D7, 10 cells): QC must exclude only D7
# (lower bound = median - 3 MAD on the log scale). The draw must be reproducible with the same
# seed, independent of the order of the cells, and never pick the excluded donor.
sizes <- c(D1 = 900, D2 = 1000, D3 = 1100, D4 = 1000, D5 = 1200, D6 = 950, D7 = 10)
toy <- do.call(rbind, lapply(names(sizes), function(d)
  data.frame(donor_id = d, cell_type = rep(c("a", "b"), length.out = sizes[[d]]))))
sp <- assign_split(toy, seed = 1)
stopifnot(
  "toy: only D7 should be excluded" = identical(sp$donor_id[sp$split == "excluded"], "D7"),
  "toy: 1 held-out and 1 E′" = sum(sp$split == "heldout") == 1 && sum(sp$split == "eprime") == 1,
  "toy: training = 6 - 2 donors" = sum(sp$split == "training") == 4,
  "toy: n_types_ge50 of D1 should be 2" = sp$n_types_ge50[sp$donor_id == "D1"] == 2,
  "toy: same seed gives the same split" = identical(assign_split(toy, seed = 1), sp),
  "toy: same input in another order gives the same split" =
    identical(assign_split(toy[nrow(toy):1, ], seed = 1), sp)
)
# Across many seeds the excluded donor is never drawn and every kept donor can be drawn.
drawn <- unlist(lapply(1:200, function(k) {
  x <- assign_split(toy, seed = k); x$donor_id[x$split %in% c("heldout", "eprime")] }))
stopifnot("toy: excluded donor drawn" = !"D7" %in% drawn,
          "toy: some kept donor is never drawn" = setequal(unique(drawn), paste0("D", 1:6)))
cat("Part 1 (toy data): PASS\n")

# Part 2: real outputs of 02_split.R
data_dir <- Sys.getenv("TESI_DATA")
split_dir <- file.path(data_dir, "derived", "split")
files <- list.files(split_dir, pattern = "_cells\\.csv\\.gz$", full.names = TRUE)
if (length(files) == 0) {
  cat("Part 2: no outputs of 02_split.R yet, skipped\n")
} else {
  for (f in files) {
    cells <- read.csv(f, colClasses = c(donor_id = "character", split = "character"))
    n_splits <- tapply(cells$split, cells$donor_id, function(s) length(unique(s)))
    stopifnot("real: a donor_id appears in two splits" = all(n_splits == 1),
              "real: a cell appears twice" = !anyDuplicated(cells$row),
              "real: cells of excluded donors in the cell table" = !any(cells$split == "excluded"),
              "real: not exactly 1 held-out and 1 E′ donor" =
                length(unique(cells$donor_id[cells$split == "heldout"])) == 1 &&
                length(unique(cells$donor_id[cells$split == "eprime"])) == 1)
    cat(sprintf("Part 2 (%s): PASS, %d donors, %d cells\n", basename(f),
                length(n_splits), nrow(cells)))
  }
}
