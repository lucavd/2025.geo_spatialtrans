# tools/S2.7_identity.R - S2.7, G-2: uguaglianza delle impronte (.md5) fra divisioni dello stesso livello e dello stesso seme.
# Uso (root del repo): Rscript --vanilla tools/S2.7_identity.R <dir> <out.csv>
a <- commandArgs(TRUE); D <- a[1]; OUT <- a[2]
f <- list.files(D, pattern = "\\.md5$", full.names = TRUE)
tag <- sub("\\.md5$", "", basename(f)); grp <- sub("_(one|p[0-9]+w[0-9]+)$", "", tag); des <- sub("^.*_", "", tag)
rows <- list()
for (g in unique(grp)) {
  k <- which(grp == g); ref <- k[des[k] == "one"]; stopifnot(length(ref) == 1L)
  H0 <- jsonlite::fromJSON(f[ref])
  for (i in setdiff(k, ref)) {
    H <- jsonlite::fromJSON(f[i])
    for (nm in names(H0)) rows[[length(rows) + 1L]] <- data.frame(layer = g, design = des[i], component = nm, md5_one = H0[[nm]], md5_design = H[[nm]], identical = identical(H0[[nm]], H[[nm]]))
  }
}
R <- do.call(rbind, rows); utils::write.csv(R, OUT, row.names = FALSE)
cat(sprintf("G-2: %d confronti, %d identici, %d diversi\n", nrow(R), sum(R$identical), sum(!R$identical)))
print(table(R$layer, R$design))
if (any(!R$identical)) { print(R[!R$identical, ]); quit(status = 1L) }
