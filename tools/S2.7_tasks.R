# tools/S2.7_tasks.R - S2.7: liste dei compiti con i semi dichiarati (docs/steps/S2.7.md, «Attese pre-registrate»). Uso (root del repo):
#   Rscript --vanilla tools/S2.7_tasks.R <G2|P1|SEC> <tasks.json>
a <- commandArgs(TRUE); Q <- a[1]; OUT <- a[2]
M <- "/mnt/micron/geo_spatialtrans"; S7 <- file.path(M, "S2.7"); RUN <- "Rscript --vanilla tools/S2.7_run.R"
draws <- function(s, n) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(s); sample.int(1e8, n) }
ENV <- "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1; "
tk <- list()
add <- function(name, cmd, cost, mem, done, log) tk[[length(tk) + 1L]] <<- list(name = name, cmd = paste0(ENV, cmd), cost = cost, mem = mem, done = done, log = log)
if (Q == "G2") {
  SQ <- c(A1 = "pc_A1_2000um_px2_s28694103", A2 = "pc_A2_2000um_px2_s79278871", A3 = "pc_A3_2000um_px2_s80792371", A4 = "pc_A4_2000um_px2_s81709255",
          A5 = "pc_A5_2000um_px2_s21892866", A6 = "pc_A6_2000um_px2_s84724683")
  s <- draws(20261218L, 7L); D <- file.path(S7, "G2"); dir.create(D, recursive = TRUE, showWarnings = FALSE)
  MEM1 <- c(A1 = 22, A2 = 32, A3 = 18, A4 = 16, A5 = 24, A6 = 8)
  for (i in seq_along(SQ)) {
    A <- names(SQ)[i]; L <- file.path(M, "S2.6/layers", paste0(SQ[[A]], ".rds"))
    for (d in list(c("one", "Inf", "1"), c("p504w4", "504", "4"), c("p1008w4", "1008", "4"))) {
      tg <- paste0("G2_", A, "_", d[1])
      add(tg, sprintf("%s new %s %d %s %s %s %s 1 full", RUN, L, s[i], d[2], d[3], D, tg), as.integer(d[3]), if (d[1] == "one") MEM1[[A]] else 4 * 8,
          file.path(D, paste0(tg, ".json")), file.path(D, paste0(tg, ".log")))
    }
  }
  L <- file.path(M, "S2.6/sections/sec_A6_s29417400.rds")
  for (d in list(c("one", "Inf", "1", 30), c("p1008w8", "1008", "8", 40), c("p1008w3", "1008", "3", 25), c("p504w8", "504", "8", 40))) {
    tg <- paste0("G2_secA6_", d[1])
    add(tg, sprintf("%s new %s %d %s %s %s %s 0 full", RUN, L, s[7], d[2], d[3], D, tg), as.integer(d[3]), as.numeric(d[4]), file.path(D, paste0(tg, ".json")), file.path(D, paste0(tg, ".log")))
  }
} else if (Q == "P1") {
  AR <- c(A1 = "A1_r1", A2 = "A2_r1", A3 = "A3_r1", A4 = "A4_f1", A5 = "A5_r1", A6 = "A6_r1"); NS <- c(A1 = 9, A2 = 18, A3 = 9, A4 = 9, A5 = 9, A6 = 27)
  S <- array(draws(20261219L, 27L * 2L * 6L), c(27L, 2L, 6L))
  Do <- file.path(S7, "p1", "old"); Dn <- file.path(S7, "p1", "new"); dir.create(Do, recursive = TRUE, showWarnings = FALSE); dir.create(Dn, recursive = TRUE, showWarnings = FALSE)
  # motore di S2.5: k = 0, 1 dal pannello di S2.5 (variante main_r2: fattore di profondita' predefinito), collegati in sola lettura
  for (A in names(AR)) for (k in 0:1) for (ext in c(".json", ".bin", ".g2.bin", ".u8.bin", ".g8.bin")) {
    src <- file.path(M, "S2.5/panel/main_r2", sprintf("%s_k%d%s", AR[[A]], k, ext)); dst <- file.path(Do, basename(src))
    if (!file.exists(dst)) file.symlink(src, dst)
  }
  for (i in seq_along(AR)) {
    A <- names(AR)[i]; L <- file.path(M, "S1.6/layers", paste0(AR[[A]], "_U6.rds"))
    for (k in seq_len(NS[[A]]) - 1L) {
      add(sprintf("P1_new_%s_k%d", A, k), sprintf("%s new %s %d 504 1 %s %s 1 full %d new", RUN, L, S[k + 1L, 2L, i], Dn, AR[[A]], k), 1L, 10,
          file.path(Dn, sprintf("%s_k%d.json", AR[[A]], k)), file.path(Dn, sprintf("%s_k%d.log", AR[[A]], k)))
      if (k >= 2L) add(sprintf("P1_old_%s_k%d", A, k), sprintf("%s old %s %d %d %s", RUN, AR[[A]], S[k + 1L, 1L, i], k, Do), 1L, 10,
                       file.path(Do, sprintf("%s_k%d.json", AR[[A]], k)), file.path(Do, sprintf("%s_k%d.log", AR[[A]], k)))
    }
  }
} else if (Q == "SEC") {
  SC <- c(A1 = "sec_A1_s84078822.rds", A2 = "sec_A2_s17624302.rds", A4 = "sec_A4_s34235060_orig.rds", A5 = "sec_A5_s14567221.rds", A6 = "sec_A6_s29417400.rds")
  s <- draws(20261220L, 5L); D <- file.path(S7, "sections"); dir.create(D, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(SC)) {
    A <- names(SC)[i]; tg <- paste0("sec_", A, "_p1008w8")
    add(tg, sprintf("%s new %s %d 1008 8 %s %s 0 light", RUN, file.path(M, "S2.6/sections", SC[[A]]), s[i], D, tg), 8L, 80,
        file.path(D, paste0(tg, ".json")), file.path(D, paste0(tg, ".log")))
  }
}
writeLines(jsonlite::toJSON(tk, auto_unbox = TRUE, pretty = TRUE), OUT)
cat(length(tk), "compiti\n")
