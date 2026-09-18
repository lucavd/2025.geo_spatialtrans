# S0 (2026-09-18): il progetto ha un DESCRIPTION, quindi renv/activate.R collocherebbe la libreria
# in ~/.cache/R/renv/library/<progetto>-<hash>. Forziamo la libreria di progetto in renv/library/
# (la stessa costruita da tools/setup_r_library.R e verificata da R/testing/test_S0_library.R).
Sys.setenv(RENV_PATHS_LIBRARY = file.path(getwd(), "renv/library"))
source("renv/activate.R")
