
M <- "/mnt/micron/geo_spatialtrans"; S7 <- file.path(M, "S2.7"); RUN <- "Rscript --vanilla tools/S2.7_run.R"; ENV <- "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; "
draws <- function(s, n) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(s); sample.int(1e8, n) }
S <- array(draws(20261219L, 27L * 2L * 6L), c(27L, 2L, 6L))
P <- file.path(S7, "p1ph"); for (d in c("old", "new")) dir.create(file.path(P, d), recursive = TRUE, showWarnings = FALSE)
for (d in c("old", "new")) for (f in list.files(file.path(S7, "p1", d), pattern = "^A5_r1_k[0-8]\\.(json|bin|g2\\.bin|u8\\.bin|g8\\.bin)$", full.names = TRUE)) {
  dst <- file.path(P, d, basename(f)); if (!file.exists(dst)) file.symlink(normalizePath(f), dst) }
tk <- list(); add <- function(n, cmd, cost, mem, done) tk[[length(tk) + 1L]] <<- list(name = n, cmd = paste0(ENV, cmd), cost = cost, mem = mem, done = done, log = sub("json$", "log", done))
L <- file.path(M, "S1.6/layers/A5_r1_U6.rds")
for (k in 9:26) {
  add(sprintf("PH_new_A5_k%d", k), sprintf("%s new %s %d 504 1 %s A5_r1 1 full %d new", RUN, L, S[k + 1L, 2L, 5L], file.path(P, "new"), k), 1L, 10, file.path(P, "new", sprintf("A5_r1_k%d.json", k)))
  add(sprintf("PH_old_A5_k%d", k), sprintf("%s old A5_r1 %d %d %s", RUN, S[k + 1L, 1L, 5L], k, file.path(P, "old")), 1L, 10, file.path(P, "old", sprintf("A5_r1_k%d.json", k)))
}
s <- draws(20261220L, 5L); D <- file.path(S7, "sections_ph")
for (x in list(c("A2", "sec_A2_s17624302.rds", 2), c("A5", "sec_A5_s14567221.rds", 4))) {
  tg <- paste0("sec_", x[1], "_p504w8"); dir.create(D, showWarnings = FALSE)
  add(tg, sprintf("%s new %s %d 504 8 %s %s 0 light", RUN, file.path(M, "S2.6/sections", x[2]), s[as.integer(x[3])], D, tg), 8L, 80, file.path(D, paste0(tg, ".json")))
}
writeLines(jsonlite::toJSON(tk, auto_unbox = TRUE, pretty = TRUE), file.path(S7, "tasks_PH.json")); cat(length(tk), "compiti\n")
