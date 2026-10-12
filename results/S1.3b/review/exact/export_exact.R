
# export esatto delle coppie P1 (revisore exact). Uso: root del repo.
REAL_ARGS <- commandArgs(TRUE)
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
suppressMessages(source("tools/S1.3_roi.R")); rm(commandArgs); source("tools/S1.3b_overlay.R")
src <- readLines("tools/S1.3b_run.R"); eval(parse(text = src[grep("^regen <- function", src):(grep("^reg_features <- function", src) - 1)]))
library(parallel)
OUT <- "/mnt/micron/geo_spatialtrans/S1.3b/review_exact"; dir.create(OUT, showWarnings = FALSE)
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)
P1 <- read.csv("results/S1.3b/S1.3b_posthoc_P1_pairs_v0.csv")
cases <- rbind(P1[, c("archetype","roi_id","model","rep","i","j")],
               data.frame(archetype="A5", roi_id="r5", model="RSA", rep=15L, i=NA, j=NA))
geom_rows <- function(g, case, role) {
  cl <- class(g)[2]; parts <- if (cl == "MULTIPOLYGON") unclass(g) else list(unclass(g))
  do.call(rbind, lapply(seq_along(parts), function(p) do.call(rbind, lapply(seq_along(parts[[p]]), function(r) {
    m <- parts[[p]][[r]]
    data.frame(case = case, role = role, type = cl, part = p, ring = r, k = seq_len(nrow(m)),
               x = sprintf("%a", m[, 1]), y = sprintf("%a", m[, 2]), xd = sprintf("%.17g", m[, 1]), yd = sprintf("%.17g", m[, 2]))
  }))))
}
one <- function(q) {
  cs <- cases[q, ]; G <- regen(roi_i(cs$archetype, cs$roi_id), cs$model, cs$rep); t <- G$terr; reg <- G$reg
  case <- sprintf("%s_%s_%s_%02d", cs$archetype, cs$roi_id, cs$model, cs$rep)
  if (!is.na(cs$i)) {
    i <- cs$i; j <- cs$j
    gx <- st_intersection(t[i], t[j])
    rows <- rbind(geom_rows(t[[i]], case, "i"), geom_rows(t[[j]], case, "j"), geom_rows(reg[[1]], case, "region"))
    if (length(gx)) rows <- rbind(rows, geom_rows(st_cast(st_collection_extract(gx, "POLYGON"), "MULTIPOLYGON")[[1]], case, "sf_inter"))
    meta <- data.frame(case = case, i = i, j = j, relate = st_relate(t[i], t[j])[1, 1],
      area_inter_sf = sprintf("%.17g", sum(as.numeric(st_area(gx)))), area_i_sf = sprintf("%.17g", as.numeric(st_area(t[i]))),
      area_j_sf = sprintf("%.17g", as.numeric(st_area(t[j]))), inter_type = paste(class(gx[[1]]), collapse = "/"),
      wkb_i = paste(as.character(st_as_binary(t[i])[[1]]), collapse = ""), wkb_j = paste(as.character(st_as_binary(t[j])[[1]]), collapse = ""),
      valid_i = st_is_valid(t[i]), valid_j = st_is_valid(t[j]), area_reg = sprintf("%.17g", G$w$area_um2),
      overlaps = st_overlaps(t[i], t[j], sparse = FALSE)[1, 1], touches = st_touches(t[i], t[j], sparse = FALSE)[1, 1])
    list(rows = rows, meta = meta)
  } else {
    rows <- rbind(geom_rows(reg[[1]], case, "region"), do.call(rbind, lapply(seq_along(t), function(k) geom_rows(t[[k]], case, sprintf("c%05d", k)))))
    a <- as.numeric(st_area(t)); u <- st_union(t)
    meta <- data.frame(case = case, n = length(t), sum_area = sprintf("%.17g", sum(a)), area_union = sprintf("%.17g", as.numeric(st_area(u))),
      area_reg_w = sprintf("%.17g", G$w$area_um2), area_reg_sf = sprintf("%.17g", as.numeric(st_area(reg))),
      symd = sprintf("%.17g", as.numeric(st_area(st_sym_difference(u, reg)))), n_invalid = sum(!st_is_valid(t)),
      union_type = class(u[[1]])[2], union_valid = st_is_valid(u))
    list(rows = rows, meta = meta, union = u)
  }
}
res <- mclapply(seq_len(nrow(cases)), function(q) tryCatch(one(q), error = function(e) conditionMessage(e)), mc.cores = 11, mc.preschedule = FALSE)
for (q in seq_along(res)) if (is.character(res[[q]])) cat("ERR", q, res[[q]], "\n")
pr <- res[seq_len(nrow(P1))]; un <- res[[nrow(cases)]]
write.csv(do.call(rbind, lapply(pr, `[[`, "rows")), file.path(OUT, "P1_pairs_coords_hex.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(pr, `[[`, "meta")), file.path(OUT, "P1_pairs_sf_meta.csv"), row.names = FALSE)
write.csv(un$rows, file.path(OUT, "A5_r5_RSA_15_coords_hex.csv"), row.names = FALSE)
write.csv(un$meta, file.path(OUT, "A5_r5_RSA_15_sf_meta.csv"), row.names = FALSE)
writeLines(paste(names(sf_extSoftVersion()), sf_extSoftVersion()), file.path(OUT, "sf_versions.txt"))
print(do.call(rbind, lapply(pr, `[[`, "meta"))[, c("case","relate","area_inter_sf","area_i_sf","inter_type","overlaps","valid_i","valid_j")]); print(un$meta)
