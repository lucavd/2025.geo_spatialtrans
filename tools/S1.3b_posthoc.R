# tools/S1.3b_posthoc.R — S1.3b: analisi POST HOC dichiarate (2026-10-10, dopo il run ufficiale, prima del report).
# P1: le coppie V0 oltre 1e-6·(Pi+Pj) (10 invece della 1 prevista) sono errori di st_intersection o sovrapposizioni reali?
#     2000 punti nell'intersezione dichiarata, appartenenza con il predicato st_intersects; D-0 sulla coppia; area V2 della coppia.
# P2: geometrie invalide di V1 e V3 (C-4): quali celle, motivo (st_is_valid reason), stato in V0.
# P3: CP-2, eccezioni GEOS di V0 vs input perturbati invalidi (stesso nullo e δ).
REAL_ARGS <- commandArgs(TRUE)
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
suppressMessages(source("tools/S1.3_roi.R")); rm(commandArgs); source("tools/S1.3b_overlay.R")
src <- readLines("tools/S1.3b_run.R"); eval(parse(text = src[grep("^regen <- function", src):(grep("^reg_features <- function", src) - 1)]))
OUTB <- if (length(REAL_ARGS)) REAL_ARGS[1] else "/mnt/micron/geo_spatialtrans/S1.3b"; RES <- "results/S1.3b"
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)
nf <- list.files(file.path(OUTB, "null"), pattern = "\\.rds$", full.names = TRUE)
hits <- list(); inval <- list()
for (f in nf) { o <- readRDS(f)
  pr <- o$res$V0$pairs; k <- which(pr$area > 1e-6 * pr$psum)
  if (length(k)) hits[[length(hits) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep, i = pr$i[k], j = pr$j[k], area_v0 = pr$area[k],
                                                        union_fail_v0 = o$res$V0$ov / o$A > 1e-9 | o$res$V0$symd / o$A > 1e-9)
  for (v in c("V1", "V3")) if (o$res[[v]]$n_invalid > 0) inval[[length(inval) + 1]] <- data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = o$rep, variant = v, n_invalid = o$res[[v]]$n_invalid)
}
hits <- do.call(rbind, hits); inval <- do.call(rbind, inval)
cases <- unique(rbind(hits[, c("archetype", "roi_id", "model", "rep")], inval[, c("archetype", "roi_id", "model", "rep")]))
P1 <- list(); P2 <- list()
for (q in seq_len(nrow(cases))) { cs <- cases[q, ]; G <- regen(roi_i(cs$archetype, cs$roi_id), cs$model, cs$rep); t <- G$terr
  h <- hits[hits$archetype == cs$archetype & hits$roi_id == cs$roi_id & hits$model == cs$model & hits$rep == cs$rep, , drop = FALSE]
  if (nrow(h)) { V2 <- s3b_variant(t, G$reg, "V2")
    for (r in seq_len(nrow(h))) { i <- h$i[r]; j <- h$j[r]
      g <- st_intersection(t[i], t[j]); set.seed(1); s <- st_sample(st_collection_extract(g, "POLYGON"), 2000)
      hi <- lengths(st_intersects(s, t[i])) > 0; hj <- lengths(st_intersects(s, t[j])) > 0
      d0 <- s3b_node(list(t[[i]], t[[j]]), 1e-6)
      P1[[length(P1) + 1]] <- data.frame(h[r, ], n_pts = length(s), in_both = sum(hi & hj), in_i_only = sum(hi & !hj), in_j_only = sum(hj & !hi),
        relate = st_relate(t[i], t[j])[1, 1], d0_diff = d0$info[["n_pairs_cross_diff"]], d0_ident = d0$info[["n_pairs_cross_identical"]], d0_tj = d0$info[["n_tjunction"]],
        d0_dmax = d0$info[["d_max"]], clipped_i = G$s$df$clipped[i], clipped_j = G$s$df$clipped[j],
        area_i = as.numeric(st_area(t[i])), area_j = as.numeric(st_area(t[j])),
        area_v2 = geos_area(geos_intersection_prec(V2$terr[i], V2$terr[j], 1e-6))) } }
  iv <- inval[inval$archetype == cs$archetype & inval$roi_id == cs$roi_id & inval$model == cs$model & inval$rep == cs$rep, , drop = FALSE]
  for (v in unique(iv$variant)) { V <- s3b_variant(t, G$reg, v); bad <- which(!st_is_valid(V$terr))
    for (b in bad) P2[[length(P2) + 1]] <- data.frame(cs, variant = v, cell = b, reason = st_is_valid(V$terr[b], reason = TRUE), valid_v0 = st_is_valid(t[b]),
      area_v0 = as.numeric(st_area(t[b])), n_vert_v0 = nrow(st_coordinates(t[b])), n_vert_v = nrow(st_coordinates(V$terr[b])), clipped = G$s$df$clipped[b],
      min_edge_v0 = { m <- st_coordinates(t[b])[, 1:2]; min(sqrt(rowSums(diff(m)^2))) }) }
}
P1 <- do.call(rbind, P1); P2 <- do.call(rbind, P2)
write.csv(P1, file.path(RES, "S1.3b_posthoc_P1_pairs_v0.csv"), row.names = FALSE); write.csv(P2, file.path(RES, "S1.3b_posthoc_P2_invalid.csv"), row.names = FALSE)
ct <- read.csv(file.path(RES, "S1.3b_cp2.csv")); c0 <- ct[ct$variant == "V0", ]
P3 <- data.frame(delta = c0$delta, exception = !is.na(c0$error), invalid_input = c0$n_invalid_input > 0, silent_fail = c0$fail & is.na(c0$error))
P3t <- aggregate(cbind(n = 1, exception, invalid_input, both = exception & invalid_input, silent_fail) ~ delta, P3, sum)
write.csv(P3t, file.path(RES, "S1.3b_posthoc_P3_cp2_exceptions.csv"), row.names = FALSE)
print(P1[, c("archetype", "roi_id", "model", "rep", "i", "j", "area_v0", "union_fail_v0", "in_both", "in_i_only", "in_j_only", "d0_diff", "d0_tj", "d0_dmax", "clipped_i", "clipped_j", "area_v2")], row.names = FALSE)
print(P2, row.names = FALSE); print(P3t, row.names = FALSE)
