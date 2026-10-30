#!/usr/bin/env Rscript
# tools/restore_library.R - S1.6b: ricostruisce da zero la libreria R del progetto dal lockfile (renv.lock), anche in un clone nuovo.
# Uso (root del repo): Rscript --vanilla tools/restore_library.R [<libreria>]
#   <libreria> predefinita: renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu (quella che usano codice, test e .Rprofile).
# Passi: (1) se renv manca, lo installa nella libreria vuota dai binari Posit PPM (Ubuntu noble); (2) renv::restore() alle versioni del
# lockfile; (3) controllo versione per versione (come tools/S1.6_env_check.R restore, nota D-1 di S1.6: i pacchetti raccomandati gia' nella
# libreria di sistema di R alla versione registrata non si reinstallano). Esce con stato 1 se un pacchetto manca o ha un'altra versione.
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
stopifnot(file.exists("renv.lock"), file.exists("DESCRIPTION"))
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
        # senza questa intestazione PPM serve i sorgenti e tutto si compila (nota del server, 2026-08-05)
        HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)),
        Ncpus = 16, install.packages.compile.from.source = "never", renv.config.pak.enabled = FALSE)
.libPaths(lib)
t0 <- Sys.time()
if (!requireNamespace("renv", lib.loc = lib, quietly = TRUE)) install.packages("renv", lib = lib)
renv::restore(project = ".", lockfile = "renv.lock", library = lib, prompt = FALSE, clean = FALSE, repos = getOption("repos"))
lk <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)$Packages
ip <- installed.packages(lib.loc = lib); is_ <- installed.packages(lib.loc = .Library)
want <- vapply(lk, function(p) p$Version, ""); have <- ip[match(names(want), rownames(ip)), "Version"]
hsys <- is_[match(names(want), rownames(is_)), "Version"]; psys <- is_[match(names(want), rownames(is_)), "Priority"]
insys <- is.na(have) & !is.na(hsys) & hsys == want & !is.na(psys) & psys == "recommended"
bad <- names(want)[!(insys | (!is.na(have) & have == want))]
cat(sprintf("%s restore: %d pacchetti nel lockfile; %d in %s alla versione registrata; %d raccomandati nella libreria di sistema di R alla versione registrata; %d diversi o mancanti; %.1f min\n",
            if (!length(bad)) "PASS" else "FAIL", length(want), sum(!is.na(have) & have == want), lib, sum(insys), length(bad),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (length(bad)) print(data.frame(pkg = bad, lock = want[bad], installed = have[match(bad, names(want))]))
quit(status = if (length(bad)) 1 else 0)
