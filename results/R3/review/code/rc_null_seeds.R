.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$i <- seq_len(nrow(rois))
f <- list.files("/mnt/micron/geo_spatialtrans/R3/null", "rds$", full.names=TRUE)
d <- do.call(rbind, lapply(f, function(p){o<-readRDS(p); data.frame(A=o$archetype, roi=o$roi_id, model=o$model, rep=o$rep, seed=o$seed, nf=o$n_failed, cv=o$summary[["cv"]], cvl=o$summary[["cv_loc"]], n=o$summary[["n"]])}))
d$i <- rois$i[match(paste(d$A,d$roi), paste(rois$archetype, rois$roi_id))]
d$exp <- 20261008 + 1000*d$i + 100*ifelse(d$model=="CSR",1,2) + d$rep
cat("n files", nrow(d), " seed==formula", sum(d$seed==d$exp), " unique seeds", length(unique(d$seed)), " unique cv", length(unique(round(d$cv,12))), "\n")
cat("RSA n_failed total", sum(d$nf[d$model=="RSA"]), " max", max(d$nf), "\n")
rr <- read.csv("results/R3/R3_roi_summary.csv")
P <- c(A1="cellpose_rgb",A2="cellpose_rgb",A3="cellpose_rgb",A4="spaceranger",A5="cellpose_rgb",A6="spaceranger")
rp <- rr[rr$method==P[rr$archetype],]
d$nobs <- rp$n[match(paste(d$A,d$roi), paste(rp$archetype, rp$roi_id))]
cat("null n == n_obs (CSR):", sum(d$n[d$model=="CSR"]==d$nobs[d$model=="CSR"]), "/", sum(d$model=="CSR"), "; RSA n-nobs range:", range(d$n[d$model=="RSA"]-d$nobs[d$model=="RSA"]), "\n")
x <- aggregate(cbind(dn=d$n-d$nobs) ~ A+model, d, function(v) round(mean(v),1)); print(x)
