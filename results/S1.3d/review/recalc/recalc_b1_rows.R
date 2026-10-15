.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

ref <- read.csv("results/S1.3/S1.3_verdicts.csv", stringsAsFactors = FALSE)
for (cd in c("V0","F1","F2","F3","F4","F5")) { x <- read.csv(file.path("/mnt/micron/geo_spatialtrans/S1.3d/b1", sprintf("verdicts_%s.csv", cd)), stringsAsFactors = FALSE)
  cols <- c("value","outcome","confirmed"); d <- which(!mapply(function(a,b) identical(a,b), split(x[, cols], seq_len(nrow(x))), split(ref[, cols], seq_len(nrow(ref)))))
  cat(cd, ": righe diverse", length(d), " sezioni", paste(unique(ref$section[d]), collapse=","), "\n")
  if (length(d)) print(data.frame(sec = ref$section[d], id = ref$id[d], ref = substr(ref$value[d],1,40), cand = substr(x$value[d],1,40), out_ref = ref$outcome[d], out_cand = x$outcome[d])) }
