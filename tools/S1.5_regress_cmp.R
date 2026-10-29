# tools/S1.5_regress_cmp.R - S1.5 (G-1): gli esiti riscritti dai test degli step passati sono identici a quelli committati?
# Confronta ogni CSV con la versione di HEAD, escluse le colonne di tempo e memoria (nomi con elapsed, time, sec, rss, mem).
# Uso: Rscript --vanilla tools/S1.5_regress_cmp.R <file.csv> [...]; stampa PASS/FAIL per file; exit 1 se qualcuno differisce.
# S1.6 (BL-165): un argomento <file versionato>=<file nuovo> confronta HEAD:<file versionato> con <file nuovo> (uscite dei test nella cartella dei log).
args <- commandArgs(TRUE); bad <- 0L
for (a in args) {
  pp <- strsplit(a, "=", fixed = TRUE)[[1]]; f <- pp[1]; nf <- if (length(pp) > 1) pp[2] else pp[1]
  old <- tryCatch(read.csv(text = system2("git", c("show", paste0("HEAD:", f)), stdout = TRUE), stringsAsFactors = FALSE), error = function(e) NULL)
  new <- tryCatch(read.csv(nf, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(old) || is.null(new)) { cat("FAIL regressione", f, "non leggibile\n"); bad <- bad + 1L; next }
  drop <- grepl("elapsed|time|sec|rss|mem|_s$", names(old), ignore.case = TRUE)
  k <- names(old)[!drop]
  # righe che misurano un tempo (testo con "tempo" o "time"): si confronta solo l'esito, non il valore misurato
  trow <- apply(old[, vapply(old, is.character, TRUE), drop = FALSE], 1, function(z) any(grepl("tempo|time", z, ignore.case = TRUE)))
  oo <- old; nn <- new
  if (any(trow) && "status" %in% k && identical(dim(old), dim(new))) { for (cc in setdiff(k, "status")) { oo[trow, cc] <- NA; nn[trow, cc] <- NA } }
  same <- identical(names(old), names(new)) && nrow(old) == nrow(new) && isTRUE(all.equal(oo[, k, drop = FALSE], nn[, k, drop = FALSE], tolerance = 1e-12, check.attributes = FALSE))
  st <- if ("status" %in% k) paste(names(table(new$status)), table(new$status), collapse = " ") else ""
  cat(sprintf("%s regressione %s identico a HEAD (colonne confrontate %d, escluse %d; righe di tempo confrontate solo sull'esito %d) %s\n", if (same) "PASS" else "FAIL", f, length(k), sum(drop), sum(trow), st))
  if (!same) bad <- bad + 1L
}
quit(status = if (bad) 1 else 0)
