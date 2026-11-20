# tools/riesame_xen_params.R - Riesame Parte B, bracci XC/XC2: legge dal pacchetto (senza modificarlo) i parametri del preset per A2 e A3
# (transcript_params, capture_params, stripe_ceiling con q = 1e-5) e li scrive in JSON. Con 'stripes' genera i fattori di striscia con
# simulate_stripes() sulla griglia vera di ogni area (ingresso: <dir>/stripes_in_<area>.parquet con row, col, in_regions; semi in
# <dir>/stripes_seeds.csv). Uso: Rscript --vanilla tools/riesame_xen_params.R <out_json> [stripes <dir>]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(jsonlite) })
for (f in sort(list.files("R", pattern = "^05b.*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE)
out <- list()
for (A in c("A2", "A3")) {
  tp <- transcript_params(A); cp <- capture_params(A)
  out[[A]] <- list(transcript = tp[c("r", "sigma_in", "sigma_e", "D", "w", "ell", "field_s2", "field_L")],
                   capture = cp[c("parity_row", "parity_col", "sd_alpha", "sd_beta", "rho_alpha", "rho_beta", "b", "w2", "ell2", "section")],
                   slow_null = is.null(cp$slow), parity_field_null = is.null(cp$parity_field),
                   Mstar = stripe_ceiling(cp, 1e-5))
}
write_json(out, a[1], auto_unbox = TRUE, digits = NA, pretty = TRUE)
if (length(a) >= 3 && a[2] == "stripes") {
  suppressPackageStartupMessages(library(arrow))
  d <- a[3]; sd <- read.csv(file.path(d, "stripes_seeds.csv"), stringsAsFactors = FALSE)
  for (h in seq_len(nrow(sd))) {
    fo <- file.path(d, sprintf("stripes_%s_k%d.parquet", sd$area[h], sd$k[h])); if (file.exists(fo)) next
    b <- as.data.frame(read_parquet(file.path(d, sprintf("stripes_in_%s.parquet", sd$area[h]))))
    r0 <- min(b$row) - min(b$row) %% 2L; c0 <- min(b$col) - min(b$col) %% 2L
    i <- b$row - r0; j <- b$col - c0; nr <- max(i) + 1L; nc <- max(j) + 1L
    tis <- logical(nr * nc); tis[i * nc + j + 1L] <- b$in_regions
    cp <- capture_params(sd$A[h])
    St <- simulate_stripes(list(grid = list(n_row = nr, n_col = nc, bin_um = 2), bins = list(in_tissue = tis)), cp, seed = sd$seed[h], tissue = tis)
    cc <- St$c[i * nc + j + 1L]
    write_parquet(data.frame(bin_id = b$bin_id, c = cc), paste0(fo, ".tmp")); file.rename(paste0(fo, ".tmp"), fo)
    cat(sd$area[h], sd$k[h], "seed", sd$seed[h], "grid", nr, nc, "norm", St$norm, "max", max(cc), "mean_tis", mean(St$c[tis]), "\n")
  }
}
