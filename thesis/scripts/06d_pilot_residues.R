# 06d_pilot_residues.R — numerical residues of the raw pilot estimates (ADEMP §7.2, item 1).
#
# Pre-registered check (ADEMP §7.2, 2026-11-03; script written 2026-11-16 before the lymphoid and brain
# runs): the residues measured in the colon pilot (sums off by at most 9e-16; RCTD 2 negative values of
# 1,079,892, the smallest -7.0e-7) are measured on the lymphoid and brain pilots. If any residue exceeds
# the tolerance 1e-5, the tolerance is revised by a dated amendment before any result.
#
# A value is a residue if it is negative but above -0.01; a spot sum is a residue if |sum - 1| <= 0.01.
# Larger deviations are failures (e.g. spots returned as all zeros), counted apart, never residues.
# PASS = no value in (-0.01, -1e-5) and no spot with |sum - 1| in (1e-5, 0.01]: every estimate is then
# either within the tolerance or clearly a failure, so the tolerance 1e-5 classifies the pilot.
# Reads the raw estimates saved by 06_pilot.R (props, before clean_props()).
#
# Output: derived/pilot/<run_id>_pilot_residues.csv, one row per method x dataset.
# Usage (working directory = thesis/):
#   Rscript scripts/06d_pilot_residues.R <tissue> [rule | rule_hd]
#   (2nd argument "rule" reads <tissue>_crule/, "rule_hd" <tissue>_crule_hd/, as 06b_pilot_nsim.R; without it, <tissue>/)

args <- commandArgs(trailingOnly = TRUE)
stopifnot("Usage: 06d_pilot_residues.R <tissue> [rule]" = length(args) >= 1)
tissue <- args[1]
run_id <- if (length(args) < 2) tissue else switch(args[2], rule = paste0(tissue, "_crule"),
                                                    rule_hd = paste0(tissue, "_crule_hd"), stop("2nd argument: rule or rule_hd"))
data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
pilot_dir <- file.path(data_dir, "derived", "pilot")
files <- list.files(file.path(pilot_dir, run_id), pattern = "^c\\d+_d\\d+_r\\d+\\.rds$", full.names = TRUE)
stopifnot("no pilot replicate files" = length(files) > 0)

tol <- 1e-5; fail_lim <- 0.01
acc <- list()
for (f in files) {
  x <- readRDS(f)
  for (d in names(x$props)) for (m in names(x$props[[d]])) {
    p <- x$props[[d]][[m]]
    na_row <- apply(p, 1, anyNA)
    v <- p[!na_row, , drop = FALSE]
    dev <- abs(rowSums(v) - 1)
    neg <- v[v < 0]
    k <- paste(m, d, sep = "|")
    a <- acc[[k]]
    if (is.null(a)) a <- list(method = m, dataset = d, files = 0, spots = 0, na_spots = 0, values = 0,
                              n_neg = 0, min_value = Inf, n_grey_neg = 0, n_fail_neg = 0,
                              max_dev_ok = 0, n_grey_sum = 0, n_fail_sum = 0)
    a$files <- a$files + 1
    a$spots <- a$spots + nrow(p); a$na_spots <- a$na_spots + sum(na_row); a$values <- a$values + length(v)
    a$n_neg <- a$n_neg + length(neg)
    if (length(v)) a$min_value <- min(a$min_value, min(v))
    a$n_grey_neg <- a$n_grey_neg + sum(neg < -tol & neg > -fail_lim)
    a$n_fail_neg <- a$n_fail_neg + sum(neg <= -fail_lim)
    ok <- dev <= fail_lim
    if (any(ok)) a$max_dev_ok <- max(a$max_dev_ok, max(dev[ok]))
    a$n_grey_sum <- a$n_grey_sum + sum(dev > tol & dev <= fail_lim)
    a$n_fail_sum <- a$n_fail_sum + sum(dev > fail_lim)
    acc[[k]] <- a
  }
}
res <- do.call(rbind, lapply(acc, as.data.frame))
res$min_value[!is.finite(res$min_value)] <- NA
res <- res[order(res$method, res$dataset), ]
rownames(res) <- NULL
res$pass <- res$n_grey_neg == 0 & res$n_grey_sum == 0

out <- file.path(pilot_dir, paste0(run_id, "_pilot_residues.csv"))
write.csv(res, out, row.names = FALSE)
print(res[, c("method", "dataset", "values", "n_neg", "min_value", "n_grey_neg", "n_fail_neg",
              "max_dev_ok", "n_grey_sum", "n_fail_sum", "pass")], digits = 3)
cat(sprintf("\n%s: %d files; %s values; %d negative, smallest %.2e; largest |sum - 1| on non-failed spots %.2e\n",
            run_id, length(files), format(sum(res$values), big.mark = ","), sum(res$n_neg),
            min(res$min_value, na.rm = TRUE), max(res$max_dev_ok)))
cat("Pre-registered check (no residue between 1e-5 and the failure limit 0.01):",
    if (all(res$pass)) "PASS" else "FAIL -> revise the tolerance by amendment (ADEMP §7.2)", "\n")
cat("Written:", out, "\n")
