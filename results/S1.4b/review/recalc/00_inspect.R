.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

x <- readRDS("/mnt/micron/geo_spatialtrans/S1.4b/loo/A1_cellpose_rgb_fold1.rds")
str(x, max.level=3, vec.len=4)
