# tools/S1.6_env_check.R - S1.6, cancello G-4 (BL-004, BL-008, BL-012, BL-078). Uso (root del repo, R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent):
#   Rscript --vanilla tools/S1.6_env_check.R desc            Imports di DESCRIPTION = pacchetti usati in R/*.R (renv::dependencies, senza i pacchetti base);
#                                                            elenco delle dipendenze dedotte da tutto il repository e non installate (BL-012)
#   Rscript --vanilla tools/S1.6_env_check.R snapshot        renv.lock dalla libreria del progetto (tutti i pacchetti installati, type = "all")
#   Rscript --vanilla tools/S1.6_env_check.R restore <lib>   renv::restore() del lockfile in una libreria vuota <lib> (binari Posit PPM); versioni = lockfile
#   Rscript --vanilla tools/S1.6_env_check.R test <lib>      R/testing/test_S1.6.R con la sola libreria <lib> (piu' la libreria di sistema di R)
LIB <- "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
args <- commandArgs(TRUE); st <- args[1]
ppm <- function() options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
  HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)),
  Ncpus = 16, install.packages.compile.from.source = "never", renv.config.pak.enabled = FALSE)
if (st == "desc") {
  .libPaths(LIB)
  d <- renv::dependencies(list.files("R", pattern = "\\.R$", full.names = TRUE), progress = FALSE, quiet = TRUE)
  base <- rownames(installed.packages(priority = "base"))
  used <- sort(setdiff(unique(d$Package), c(base, "geo.spatialtrans")))
  imp <- trimws(sub("\\(.*", "", strsplit(read.dcf("DESCRIPTION", fields = "Imports")[1, 1], ",")[[1]]))
  imp <- sort(imp[nzchar(imp)])
  ok <- identical(used, imp)
  cat(sprintf("%s Imports di DESCRIPTION = pacchetti usati in R/*.R (%d): %s\n", if (ok) "PASS" else "FAIL", length(used), paste(used, collapse = ", ")))
  if (!ok) cat("  solo nel codice:", setdiff(used, imp), " | solo in DESCRIPTION:", setdiff(imp, used), "\n")
  all <- renv::dependencies(".", progress = FALSE, quiet = TRUE)
  inst <- c(rownames(installed.packages(lib.loc = c(LIB, .Library))))
  miss <- sort(setdiff(unique(all$Package), c(inst, base, "geo.spatialtrans")))
  cat(sprintf("INFO dipendenze dedotte da tutto il repository e non installate (BL-012): %d\n", length(miss)))
  for (p in miss) cat("  ", p, ":", paste(unique(dirname(all$Source[all$Package == p])), collapse = ", "), "\n")
  quit(status = if (ok) 0 else 1)
}
if (st == "snapshot") {
  .libPaths(LIB)
  renv::snapshot(project = ".", library = LIB, lockfile = "renv.lock", type = "all", prompt = FALSE, force = TRUE)
  lk <- jsonlite::fromJSON("renv.lock"); cat("lockfile:", length(lk$Packages), "pacchetti\n")
}
if (st == "restore") {
  lib <- args[2]; dir.create(lib, showWarnings = FALSE, recursive = TRUE); stopifnot(length(list.files(lib)) == 0)
  ppm(); t0 <- Sys.time()
  .libPaths(c(lib, LIB))                       # renv stesso viene dalla libreria del progetto; i pacchetti vanno in <lib>
  renv::restore(project = ".", lockfile = "renv.lock", library = lib, prompt = FALSE, clean = FALSE, repos = getOption("repos"))
  lk <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)$Packages
  ip <- installed.packages(lib.loc = lib)
  want <- vapply(lk, function(p) p$Version, ""); have <- ip[match(names(want), rownames(ip)), "Version"]
  bad <- names(want)[is.na(have) | have != want]
  cat(sprintf("%s restore da zero: %d pacchetti nel lockfile, %d installati in %s alla versione registrata, %d diversi o mancanti, %.1f min\n",
              if (!length(bad)) "PASS" else "FAIL", length(want), sum(!is.na(have) & have == want), lib, length(bad), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (length(bad)) print(data.frame(pkg = bad, lock = want[bad], installed = have[match(bad, names(want))]))
  quit(status = if (length(bad)) 1 else 0)
}
if (st == "test") {
  lib <- args[2]; Sys.setenv(S16_LIB = lib)
  quit(status = system2("Rscript", c("--vanilla", "R/testing/test_S1.6.R"), env = c(paste0("S16_LIB=", lib), "R_LIBS_SITE=/nonexistent", "R_LIBS_USER=/nonexistent")))
}
