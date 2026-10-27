# 01_download.R — download of the real scRNA-seq data (PROTOCOLLO §4.2 step 7, §6 week 2)
# Downloads the three CELLxGENE h5ad files chosen in DECISIONS (row "Dataset", 2026-10-26)
# into TESI_DATA/raw, computes their md5 and writes, or checks, thesis/data_manifest.csv
# (columns file, url, bytes, md5, data_download).
#
# Pre-registered checks (written before the first run, 2026-10-27):
#   1. each file has exactly the number of bytes the CELLxGENE API reported on 2026-10-27
#      (expected_bytes below); otherwise the file is deleted and the script stops;
#   2. if the manifest already lists the file (e.g. rerun on the server, §4.2 step 9),
#      its md5 must match; a mismatch stops the script and nothing is overwritten;
#   3. the version pinned below must still be the current one on CELLxGENE; if not, the
#      script only warns (the pinned version is still downloaded, so the md5 stays comparable).
# Whole files are downloaded. Subsets (Wells thoracic LN, Pelka 3' v2 MMRp tumours,
# PsychAD healthy adults) are made in 02_split.R into data/derived/: raw files are read-only (§3.7).
#
# Run from the thesis project (working directory = thesis/):
#   source("scripts/01_download.R")       or   Rscript scripts/01_download.R
# Downloads are ~29 GB in total and are not resumed: on the server, run it inside tmux.

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
raw_dir <- file.path(data_dir, "raw")
stopifnot("Missing TESI_DATA/raw" = dir.exists(raw_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
stopifnot("Package jsonlite missing: renv::restore()" = requireNamespace("jsonlite", quietly = TRUE))
manifest_path <- "data_manifest.csv"

# Pinned dataset versions (CELLxGENE Curation API, read 2026-10-27).
# A dataset_version_id never changes content, so the PC and the server get the same bytes.
datasets <- data.frame(
  tissue = c("brain", "colon", "lymphoid"),
  collection_id = c("84ce6837-548d-4a1f-919f-0bc0d9a3952f",   # PsychAD, doi:10.1101/2024.10.31.24316513
                    "4d82bd8e-8827-410b-8fc7-685b7dd92585",   # Pelka 2021, doi:10.1016/j.cell.2021.08.003
                    "cc431242-35ea-41e1-a100-41e0dec2665b"),  # Wells 2025, doi:10.1038/s41590-025-02241-4
  dataset_title = c("Aging_Cohort", "All Cells", "Immune aging project - all cells"),
  dataset_version_id = c("13e8f1dd-962f-47b7-9cf7-c71d2b21b8a5",
                         "3d1f3641-5755-47e2-b797-10317d33fe32",
                         "589324ee-ce60-4295-843f-29a52894a522"),
  expected_bytes = c(12421671635, 3759292109, 12689317573),
  stringsAsFactors = FALSE
)
datasets$url  <- sprintf("https://datasets.cellxgene.cziscience.com/%s.h5ad", datasets$dataset_version_id)
datasets$file <- sprintf("%s_%s.h5ad", datasets$tissue, datasets$dataset_version_id)

# Check 3: is the pinned version still the current one?
check_current_version <- function(collection_id, title, version_id) {
  api <- sprintf("https://api.cellxgene.cziscience.com/curation/v1/collections/%s", collection_id)
  coll <- tryCatch(jsonlite::fromJSON(api), error = function(e) NULL)
  if (is.null(coll)) {
    warning("CELLxGENE API not reachable; version of '", title, "' not checked")
    return(invisible(NULL))
  }
  current <- coll$datasets$dataset_version_id[coll$datasets$title == title]
  if (!version_id %in% current) {
    warning("'", title, "': pinned version ", version_id, " is no longer the current one (now ",
            paste(current, collapse = ", "), "). Downloading the pinned version.")
  }
}

manifest <- if (file.exists(manifest_path)) {
  read.csv(manifest_path, colClasses = "character")
} else {
  data.frame(file = character(), url = character(), bytes = character(),
             md5 = character(), data_download = character())
}

options(timeout = 6 * 3600)   # seconds; large files on a slow connection

for (i in seq_len(nrow(datasets))) {
  d <- datasets[i, ]
  dest <- file.path(raw_dir, d$file)
  cat(sprintf("\n[%s] %s (%.2f GB)\n", d$tissue, d$dataset_title, d$expected_bytes / 1e9))
  check_current_version(d$collection_id, d$dataset_title, d$dataset_version_id)

  if (!file.exists(dest)) {
    part <- paste0(dest, ".part")
    status <- download.file(d$url, part, mode = "wb", quiet = FALSE)
    stopifnot("download.file failed" = status == 0)
    file.rename(part, dest)
  } else {
    cat("  already on disk, not downloaded again\n")
  }

  # Check 1: size
  bytes <- file.size(dest)
  if (bytes != d$expected_bytes) {
    file.remove(dest)
    stop(sprintf("%s: %.0f bytes, expected %.0f; file deleted", d$file, bytes, d$expected_bytes))
  }

  # Check 2: md5 against the manifest
  md5 <- unname(tools::md5sum(dest))
  row <- which(manifest$file == d$file)
  if (length(row) == 1) {
    if (manifest$md5[row] != md5) {
      stop(sprintf("%s: md5 %s differs from manifest %s", d$file, md5, manifest$md5[row]))
    }
    cat("  md5 matches the manifest:", md5, "\n")
  } else {
    manifest <- rbind(manifest, data.frame(
      file = d$file, url = d$url, bytes = format(bytes, scientific = FALSE),
      md5 = md5, data_download = format(Sys.Date())))
    write.csv(manifest, manifest_path, row.names = FALSE)
    cat("  added to the manifest:", md5, "\n")
  }
}

cat("\nAll files present, sizes and md5 OK. Manifest:", normalizePath(manifest_path), "\n")
