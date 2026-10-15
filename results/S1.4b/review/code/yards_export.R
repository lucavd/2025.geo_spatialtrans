
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
y <- do.call(rbind, lapply(list.files("/tmp/code_vc", pattern = "^A[24]_fold[0-9]\\.rds$", full.names = TRUE), function(f) { o <- readRDS(f)
  m <- o$m[o$m$variant == "V-ab", ]
  data.frame(file = sprintf("%s_fold%d_s1.parquet", o$A, o$k), y_nc_med = median(o$yard$ks_nc), y_nc_max = max(o$yard$ks_nc), y_an_med = median(o$yard$ks_an), y_an_max = max(o$yard$ks_an),
             y_sp_med = median(o$yard$d_sp), y_sp_max = max(o$yard$d_sp), ab_ks_nc_med10 = median(m$ks_nc), ab_ks_an_med10 = median(m$ks_an), ab_ks_an_s1 = m$ks_an[m$seed_idx == 1]) }))
write.csv(y, "/tmp/code_vc/yards_A2A4.csv", row.names = FALSE)
