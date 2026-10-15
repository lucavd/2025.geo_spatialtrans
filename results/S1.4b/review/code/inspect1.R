
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
o <- readRDS("/mnt/micron/geo_spatialtrans/S1.4/loo/A2_cellpose_rgb_fold1.rds")
str(o, max.level = 1)
cat("class seeds:", class(o$seeds), "\n"); print(head(o$seeds, 3)); print(o$sims$M2)
str(o$yard)
r <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.4/real/A2_%s_cellpose_rgb.rds", o$roi_test))
print(summary(r$cells[, c("x","y","a_t","a_n","d_gen","d_poi")])); cat("n", r$n, "\n")
b <- readRDS("/mnt/micron/geo_spatialtrans/S1.4b/loo/A2_cellpose_rgb_fold1.rds"); str(b, max.level=1); print(head(b$m)); print(head(b$checks,3))
