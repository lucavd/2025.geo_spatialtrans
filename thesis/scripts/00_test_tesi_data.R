# 00_test_tesi_data.R — test script (PROTOCOLLO §6, week 1)
# Checks that R finds the data folder through TESI_DATA (§4.2 step 5)
# and that its subfolders raw, derived and sim exist. Reads and writes no data.

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))

sub_names <- c("raw", "derived", "sim")
found <- dir.exists(file.path(data_dir, sub_names))
names(found) <- sub_names

cat("TESI_DATA:", data_dir, "\n")
for (s in sub_names) cat(sprintf("  %-8s %s\n", s, if (found[[s]]) "OK" else "MISSING"))
stopifnot("Missing subfolder(s) in TESI_DATA" = all(found))
cat("All paths OK.\n")
