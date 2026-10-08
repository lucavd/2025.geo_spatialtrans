.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); suppressPackageStartupMessages(library(arrow))
x <- as.data.frame(read_parquet("/mnt/micron/geo_spatialtrans/R3/cells_all.parquet", col_select = c("archetype","roi_id","method","interior","area")))
x <- x[x$archetype=="A5" & x$method=="cellpose_rgb" & x$interior,]
for (r in sort(unique(x$roi_id))) { a <- x$area[x$roi_id==r]; cat(r, length(a), diptest::dip.test(a)$p.value, diptest::dip.test(log(a))$p.value, diptest::dip.test(sqrt(a))$p.value, "\n") }

