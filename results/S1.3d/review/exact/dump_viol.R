# review/exact/dump_viol.R — revisore esatto S1.3d: estrae TUTTE le violazioni salvate in res/*.rds (sola lettura) in un parquet.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(arrow))
RES <- "/mnt/micron/geo_spatialtrans/S1.3d/res"; OUT <- "/mnt/micron/geo_spatialtrans/S1.3d/review_exact"
fs <- list.files(RES, pattern = "\\.rds$", full.names = TRUE)
L <- list(); C5 <- list(); nm_all <- character(0)
for (f in fs) { o <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(o)) next
  key <- if (!is.null(o$key)) o$key else if (grepl("^S_cperf", basename(f))) "S_cperf" else basename(f)
  v <- o$viol; if (is.null(v) && !is.null(o$cmp)) v <- list(cmp = o$cmp$viol)
  if (length(v)) for (nm in names(v)) { d <- v[[nm]]; if (is.null(d) || !is.data.frame(d) || !nrow(d)) next
      nm_all <- c(nm_all, nm)
      d <- as.data.frame(lapply(d, function(z) if (is.factor(z)) as.character(z) else z), stringsAsFactors = FALSE)
      d$key <- key; d$vname <- nm; d$file <- basename(f); L[[length(L) + 1]] <- d }
  if (!is.null(o$c5) && !is.null(o$c5$info)) C5[[length(C5) + 1]] <- data.frame(key = key, file = basename(f),
       n_invalid = if (is.null(o$c5$info$n_invalid)) NA else o$c5$info$n_invalid,
       invalid_ids = if (is.null(o$c5$info$invalid_ids)) NA_character_ else paste(o$c5$info$invalid_ids, collapse = ";"))
}
nm <- unique(unlist(lapply(L, names)))
V <- do.call(rbind, lapply(L, function(d) { for (k in setdiff(nm, names(d))) d[[k]] <- NA; d[, nm, drop = FALSE] }))
write_parquet(V, file.path(OUT, "viol_all.parquet"))
C5 <- do.call(rbind, C5); write.csv(C5, file.path(OUT, "c5_info.csv"), row.names = FALSE)
cat("files", length(fs), "viol rows", nrow(V), "\n"); print(table(V$vname))
cat("c5 n_invalid>0:", sum(C5$n_invalid > 0, na.rm = TRUE), " total", sum(C5$n_invalid, na.rm = TRUE), "\n")
print(names(o)); if (!is.null(o$c5)) str(o$c5$info)
