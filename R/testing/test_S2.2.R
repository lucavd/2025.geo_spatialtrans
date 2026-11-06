# R/testing/test_S2.2.R - S2.2: batteria di test di simulate_cell_profiles() e load_cell_profiles() (cancello G-1).
# Pre-registrazione: docs/steps/S2.2.md. Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S22_RSRC=<cartella con R/>] Rscript --vanilla R/testing/test_S2.2.R
# I mutanti sono copie modificate del codice (tools/S2.2_mutants.sh) passate con S22_RSRC.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf) })
SRC <- Sys.getenv("S22_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
err <- function(expr) inherits(tryCatch({ expr; NULL }, error = function(e) e), "error")
# --- profili sintetici: 200 geni, tipo «a» con 3 sottotipi (quote 0.5, 0.3, 0.2), tipo «b» con 1 ---
G <- 200L
mkP <- function(seed) { set.seed(seed); x <- rgamma(G, 0.3); x / sum(x) }
P <- rbind(mkP(1), mkP(2), mkP(3), mkP(4)); rownames(P) <- c("a|1", "a|2", "a|3", "b|1"); colnames(P) <- sprintf("g%03d", seq_len(G))
prof <- list(archetype = "T", genes = data.frame(gene_id = colnames(P), gene_name = colnames(P)),
             subtypes = data.frame(subtype = rownames(P), cell_type = c("a", "a", "a", "b"), share = c(0.5, 0.3, 0.2, 1), stringsAsFactors = FALSE),
             P = P, phi = rep(0, G), totals = data.frame(q = seq(0, 1, by = 0.001), a = seq(50, 150, length.out = 1001), b = seq(200, 300, length.out = 1001)),
             noise = "M", depth_factor = 1.5, meta = list(winner = "test"))
cells <- data.frame(cell_id = 1000L + seq_len(22000L), cell_type = c(rep("a", 20000L), rep("b", 2000L)), stringsAsFactors = FALSE)
s1 <- simulate_cell_profiles(cells, prof, seed = 11)
# (a) totali
ck("a1 somma dei conteggi = totale per ogni cellula", all(Matrix::colSums(s1$counts) == s1$cells$n_transcripts))
ck("a2 totali >= 1 e interi", all(s1$cells$n_transcripts >= 1) && all(s1$cells$n_transcripts == round(s1$cells$n_transcripts)))
mA <- mean(s1$cells$n_transcripts[s1$cells$cell_type == "a"]); seA <- sd(s1$cells$n_transcripts[s1$cells$cell_type == "a"]) / sqrt(20000)
ck("a3 totale medio = fattore x media della distribuzione (|z| < 5)", abs(mA - 1.5 * 100) / seA < 5, sprintf("media %.2f atteso 150", mA))
ck("a4 dimensioni e nomi", identical(dim(s1$counts), c(G, 22000L)) && identical(colnames(s1$counts), as.character(cells$cell_id)) &&
     identical(rownames(s1$counts), colnames(P)))
# (b) quote e profili ricostruiti
zmax_prof <- 0; zmax_share <- 0
for (s in rownames(P)) { w <- s1$cells$subtype == s; obs <- Matrix::rowSums(s1$counts[, w, drop = FALSE]); Tt <- sum(s1$cells$n_transcripts[w])
  ex <- Tt * P[s, ]; k <- ex >= 20; z <- (obs[k] - ex[k]) / sqrt(ex[k] * (1 - P[s, k])); zmax_prof <- max(zmax_prof, abs(z)) }
for (t in c("a")) { r <- prof$subtypes$cell_type == t; n_t <- sum(s1$cells$cell_type == t)
  for (s in prof$subtypes$subtype[r]) { p <- prof$subtypes$share[prof$subtypes$subtype == s]; o <- sum(s1$cells$subtype == s)
    zmax_share <- max(zmax_share, abs(o - n_t * p) / sqrt(n_t * p * (1 - p))) } }
ck("b1 profili ricostruiti dalla somma delle cellule (|z| < 5 per ogni gene con atteso >= 20)", zmax_prof < 5, sprintf("max |z| %.2f", zmax_prof))
ck("b2 quote dei sottotipi (|z| < 5)", zmax_share < 5, sprintf("max |z| %.2f", zmax_share))
ck("b3 sottotipi solo del proprio tipo", all(startsWith(s1$cells$subtype, paste0(s1$cells$cell_type, "|"))))
# (c) semi, processi, stato del generatore
set.seed(99); before <- .Random.seed
s2 <- simulate_cell_profiles(cells, prof, seed = 11)
ck("c1 stesso seme, stesso risultato", identical(s1$counts, s2$counts) && identical(s1$cells, s2$cells))
ck("c2 stato del generatore dell'utente ripristinato", identical(before, .Random.seed))
s4 <- simulate_cell_profiles(cells, prof, seed = 11, workers = 4L)
ck("c3 risultato identico con 1 e 4 processi", identical(s1$counts, s4$counts) && identical(s1$cells, s4$cells))
s3 <- simulate_cell_profiles(cells, prof, seed = 12)
ck("c4 seme diverso, risultato diverso", !identical(s1$counts, s3$counts))
# (d) errori espliciti
bad <- prof; bad$subtypes$cell_type[4] <- "c"; bad$subtypes$share[4] <- 1
ck("d1 tipo assente dai profili -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
bad <- prof; bad$P[1, 1] <- bad$P[1, 1] + 0.01
ck("d2 profilo che non somma a 1 -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
bad <- prof; bad$P[1, 1] <- -bad$P[1, 1]; bad$P[1, 2] <- bad$P[1, 2] + 2 * prof$P[1, 1]
ck("d3 valore negativo -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
bad <- prof; bad$P[1, 1] <- NA
ck("d4 NA nei profili -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
bad <- prof; bad$phi <- bad$phi[-1]
ck("d5 phi di lunghezza diversa dai geni -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
bad <- prof; bad$genes <- bad$genes[-1, ]
ck("d6 scala dei geni diversa dai profili -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
ck("d7 cell_df senza cell_type -> errore", err(simulate_cell_profiles(cells[, "cell_id", drop = FALSE], prof, seed = 1)))
ck("d8 senza seme -> errore", err(simulate_cell_profiles(cells, prof)))
ck("d9 fattore di profondita' <= 0 -> errore", err(simulate_cell_profiles(cells, prof, seed = 1, depth_factor = 0)))
bad <- prof; bad$totals$b <- NULL
ck("d10 tipo senza distribuzione dei totali -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
cdup <- cells; cdup$cell_id[2] <- cdup$cell_id[1]
ck("d11 cell_id duplicati -> errore", err(simulate_cell_profiles(cdup, prof, seed = 1)))
bad <- prof; bad$subtypes$share[1] <- 0.6
ck("d12 quote dei sottotipi che non sommano a 1 -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
bad <- prof; bad$totals$a[5] <- NA
ck("d13 NA nei quantili dei totali -> errore", err(simulate_cell_profiles(cells, bad, seed = 1)))
# (e) catalogo dei tipi dei preset = catalogo del banco
arch <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
py <- system2(".venv/bin/python", c("-c", shQuote("import sys; sys.path.insert(0, 'tools'); import S22_common as S\nfor a, t in S.PRESET_TYPES.items(): print(a + '\\t' + '|'.join(sorted(t)))")), stdout = TRUE)
pyt <- setNames(sub("^[^\t]*\t", "", py), sub("\t.*$", "", py))
okc <- vapply(names(arch), function(a) { r <- sort(unique(tissue_preset_spec(arch[[a]])$composition$cell_type)); identical(paste(r, collapse = "|"), pyt[[a]]) }, logical(1))
ck("e1 tipi dei preset (R) = PRESET_TYPES di tools/S22_common.py", all(okc), paste(names(okc)[!okc], collapse = ","))
rl <- readLines("tools/S2.2_real.py"); cl <- readLines("tools/S22_common.py")
dict_txt <- function(l) { i <- grep("^PRESET_TYPES = ", l); j <- i - 1L + grep("\\]\\}\\s*$", l[i:length(l)])[1]; paste(trimws(l[i:j]), collapse = " ") }
ck("e2 PRESET_TYPES identico in tools/S2.2_real.py e tools/S22_common.py", identical(dict_txt(rl), dict_txt(cl)))
# (f) modello G: varianza fra cellule di un gene = formula esatta (gamma sul solo gene, totale fisso)
pg <- prof; pg$phi <- rep(0, G); pg$phi[7] <- 0.5; pg$totals$a <- rep(1000, 1001); pg$P[1, 7] <- 0
pg$P[1, ] <- pg$P[1, ] * (1 - 0.05) / sum(pg$P[1, ]); pg$P[1, 7] <- 0.05
pg$subtypes$share <- c(1, 0, 0, 1)
cg <- data.frame(cell_id = seq_len(20000L), cell_type = "a", stringsAsFactors = FALSE)
sg <- simulate_cell_profiles(cg, pg, seed = 5, noise = "G", depth_factor = 1)
x7 <- as.numeric(sg$counts[7, ]); N <- 1000; p <- 0.05; k <- 1 / 0.5
lamf <- function(w) p * w / (p * w + 1 - p)
E1 <- integrate(function(w) lamf(w) * dgamma(w, k, k), 0, Inf)$value; E2 <- integrate(function(w) lamf(w)^2 * dgamma(w, k, k), 0, Inf)$value
vex <- N * (E1 - E2) + N^2 * (E2 - E1^2)
ck("f1 modello G: varianza del gene = formula esatta (entro 5 %)", abs(var(x7) / vex - 1) < 0.05, sprintf("osservata %.1f attesa %.1f", var(x7), vex))
sm <- simulate_cell_profiles(cg, pg, seed = 5, noise = "M", depth_factor = 1)
ck("f2 modello M sullo stesso caso: varianza binomiale (entro 5 %)", abs(var(as.numeric(sm$counts[7, ])) / (N * p * (1 - p)) - 1) < 0.05)
# (g) sottotipi spenti
sn <- simulate_cell_profiles(cells, prof, seed = 11, use_subtypes = FALSE)
pa <- colSums(P[1:3, ] * c(0.5, 0.3, 0.2)); wa <- sn$cells$cell_type == "a"; obs <- Matrix::rowSums(sn$counts[, wa]); ex <- sum(sn$cells$n_transcripts[wa]) * pa
k2 <- ex >= 20
ck("g1 sottotipi spenti: un solo profilo per tipo, media pesata dei sottotipi", all(sn$cells$subtype %in% c("a", "b")) &&
     max(abs((obs[k2] - ex[k2]) / sqrt(ex[k2]))) < 5)
# (h) lettura dei file nel formato di tools/S2.2_build_profiles.py
td <- tempfile("s22prof"); dir.create(td)
S0 <- round(P * 1e5); trip <- which(S0 > 0, arr.ind = TRUE)
write.table(data.frame(subtype = rownames(P)[trip[, 1]], gene_index = trip[, 2], count = S0[trip]), gzfile(file.path(td, "T_sums.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(subtype = rownames(P), cell_type = prof$subtypes$cell_type, n_cells = 100, share = prof$subtypes$share), file.path(td, "T_subtypes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(prof$genes, file.path(td, "T_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene_index = seq_len(G), phi = seq_len(G) / 1000), file.path(td, "T_phi.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(prof$totals, file.path(td, "T_totals.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines('{\n "archetype": "T",\n "winner": "exp_G",\n "source": "exp",\n "noise": "G",\n "depth_factor": 1.25,\n "atlas": null\n}', file.path(td, "T_meta.json"))
lp <- load_cell_profiles("T", dir = td)
ck("h1 lettura: profili = somme / totale, nell'ordine dei geni", max(abs(lp$P - S0 / rowSums(S0))) < 1e-12 && identical(colnames(lp$P), colnames(P)))
ck("h2 lettura: phi, modello, fattore", identical(lp$phi, seq_len(G) / 1000) && lp$noise == "G" && lp$depth_factor == 1.25)
# (i) totali per sottotipo
ps <- prof; ps$totals_subtype <- data.frame(q = prof$totals$q, `a|1` = seq(10, 20, length.out = 1001), `a|2` = seq(100, 200, length.out = 1001),
  `a|3` = seq(1000, 2000, length.out = 1001), `b|1` = seq(200, 300, length.out = 1001), check.names = FALSE)
si_ <- simulate_cell_profiles(cells, ps, seed = 21, totals_by = "subtype", depth_factor = 1)
mm <- tapply(si_$cells$n_transcripts, si_$cells$subtype, mean); ss <- tapply(si_$cells$n_transcripts, si_$cells$subtype, function(v) sd(v) / sqrt(length(v)))
exm <- c(`a|1` = 15, `a|2` = 150, `a|3` = 1500, `b|1` = 250)
ck("i1 totali per sottotipo: media di ogni sottotipo = media della sua distribuzione (|z| < 5)", all(abs(mm[names(exm)] - exm) / ss[names(exm)] < 5),
   paste(sprintf("%s %.1f", names(mm), mm), collapse = "; "))
ck("i2 totali per sottotipo: somma = totale, quote come prima", all(Matrix::colSums(si_$counts) == si_$cells$n_transcripts) &&
     abs(mean(si_$cells$subtype[si_$cells$cell_type == "a"] == "a|1") - 0.5) < 5 * sqrt(0.25 / 20000))
ck("i3 totali per sottotipo senza i quantili -> errore", err(simulate_cell_profiles(cells, prof, seed = 1, totals_by = "subtype")))
s5 <- simulate_cell_profiles(cells, ps, seed = 11)
ck("i4 senza totals_by nei metadati il predefinito e' \"type\" e il risultato non cambia se ci sono anche i quantili per sottotipo", identical(s5$counts, s1$counts))
ps2 <- ps; ps2$meta$totals_by <- "subtype"
s6 <- simulate_cell_profiles(cells, ps2, seed = 21, depth_factor = 1)
ck("i5 con totals_by = \"subtype\" nei metadati il predefinito segue i metadati (DEC-046)", identical(s6$counts, si_$counts) && s6$metadata$totals_by == "subtype")
cat(sprintf("\n%d PASS, %d FAIL\n", sum(unlist(res)), sum(!unlist(res))))
if (!all(unlist(res))) quit(status = 1)
