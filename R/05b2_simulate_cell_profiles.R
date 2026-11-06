# R/05b2_simulate_cell_profiles.R - Step 2, S2.2: trascritti per cellula (totale dal tipo, sottotipo, conteggi per gene).
# Design dello Step 2 v0.2 (docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md), sezione 4 e D2-3, D2-8, D2-9;
# fonte dei profili, modello di conteggio e fattore di profondita' decisi dai dati in S2.2 (docs/steps/S2.2.md, DEC-040).
# I profili si costruiscono con tools/S2.2_build_profiles.py (file su micron, md5 in results/S2.2/S2.2_profiles_md5.txt).
# Dipendenze: Matrix, parallel (base).

#' Carica l'insieme dei profili di un tessuto
#'
#' @param archetype `"A1"` ... `"A6"`.
#' @param dir cartella con i file di `tools/S2.2_build_profiles.py` (`<A>_sums.tsv.gz`, `_subtypes.tsv`, `_genes.tsv`, `_phi.tsv`,
#'   `_totals.tsv`, `_meta.json`; facoltativo `_totals_subtype.tsv`, quantili dei totali per sottotipo).
#' @return lista `archetype`, `genes` (data.frame `gene_id`, `gene_name`), `subtypes` (data.frame `subtype`, `cell_type`, `share`),
#'   `P` (matrice sottotipi x geni, righe che sommano a 1), `phi` (per gene), `totals` (data.frame: colonna `q` e una colonna per tipo,
#'   quantili dei trascritti per cellula), `totals_subtype` (come `totals`, una colonna per sottotipo, o NULL), `noise` (`"M"` multinomiale o `"G"` con fattore gamma per gene), `depth_factor`, `meta`.
#' @export
load_cell_profiles <- function(archetype, dir = getOption("geo.spatialtrans.profiles_dir", "/mnt/micron/geo_spatialtrans/S2.2/profiles")) {
  f <- function(s) file.path(dir, paste0(archetype, s))
  if (!file.exists(f("_meta.json"))) stop("load_cell_profiles(): profili assenti per ", archetype, " in ", dir, call. = FALSE)
  rd <- function(p) utils::read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
  genes <- rd(f("_genes.tsv")); sub <- rd(f("_subtypes.tsv")); sums <- rd(gzfile(f("_sums.tsv.gz")))
  phi <- rd(f("_phi.tsv")); tot <- rd(f("_totals.tsv"))
  tots <- if (file.exists(f("_totals_subtype.tsv"))) rd(f("_totals_subtype.tsv")) else NULL
  meta <- .s22_read_json(f("_meta.json"))
  S <- Matrix::sparseMatrix(i = match(sums$subtype, sub$subtype), j = sums$gene_index, x = as.numeric(sums$count),
                            dims = c(nrow(sub), nrow(genes)))
  rs <- Matrix::rowSums(S)
  if (any(rs <= 0)) stop("load_cell_profiles(): sottotipi senza conteggi: ", paste(sub$subtype[rs <= 0], collapse = ", "), call. = FALSE)
  P <- as.matrix(S / rs)
  dimnames(P) <- list(sub$subtype, genes$gene_id)
  ph <- numeric(nrow(genes)); ph[phi$gene_index] <- phi$phi
  list(archetype = archetype, genes = genes, subtypes = sub[, c("subtype", "cell_type", "share")], P = P, phi = ph, totals = tot,
       totals_subtype = tots, noise = meta$noise, depth_factor = as.numeric(meta$depth_factor), meta = meta)
}

.s22_read_json <- function(p) {
  # lettore minimale per i meta piatti di tools/S2.2_build_profiles.py (chiave: valore), senza dipendenze
  txt <- paste(readLines(p, warn = FALSE), collapse = "")
  txt <- gsub("^\\s*\\{|\\}\\s*$", "", txt)
  parts <- strsplit(txt, ",\\s*(?=\")", perl = TRUE)[[1]]
  out <- list()
  for (x in parts) {
    kv <- regmatches(x, regexec("^\\s*\"([^\"]+)\"\\s*:\\s*(.*)$", x))[[1]]
    if (length(kv) != 3) next
    v <- trimws(kv[3])
    out[[kv[2]]] <- if (grepl("^\"", v)) gsub("^\"|\"$", "", v) else if (v == "null") NA else suppressWarnings(as.numeric(v))
  }
  out
}

.s22_check_profiles <- function(prof, types) {
  P <- prof$P
  if (!is.matrix(P) || anyNA(P) || any(P < 0)) stop("simulate_cell_profiles(): profili con valori negativi o NA.", call. = FALSE)
  if (any(abs(rowSums(P) - 1) > 1e-9)) stop("simulate_cell_profiles(): un profilo non somma a 1.", call. = FALSE)
  if (length(prof$phi) != ncol(P) || anyNA(prof$phi) || any(prof$phi < 0))
    stop("simulate_cell_profiles(): `phi` non ha un valore >= 0 per ogni gene del profilo.", call. = FALSE)
  if (nrow(prof$genes) != ncol(P)) stop("simulate_cell_profiles(): geni e profili non hanno la stessa lunghezza.", call. = FALSE)
  if (!identical(rownames(P), prof$subtypes$subtype)) stop("simulate_cell_profiles(): sottotipi e righe dei profili non coincidono.", call. = FALSE)
  sh <- tapply(prof$subtypes$share, prof$subtypes$cell_type, sum)
  if (anyNA(prof$subtypes$share) || any(prof$subtypes$share < 0) || any(abs(sh - 1) > 1e-9))
    stop("simulate_cell_profiles(): le quote dei sottotipi di un tipo non sommano a 1.", call. = FALSE)
  tt <- prof$totals[, setdiff(names(prof$totals), "q"), drop = FALSE]
  if (anyNA(tt) || any(tt < 0)) stop("simulate_cell_profiles(): quantili dei totali con NA o negativi.", call. = FALSE)
  miss <- setdiff(types, prof$subtypes$cell_type)
  if (length(miss)) stop("simulate_cell_profiles(): tipi senza profilo: ", paste(miss, collapse = ", "), call. = FALSE)
  miss <- setdiff(types, setdiff(names(prof$totals), "q"))
  if (length(miss)) stop("simulate_cell_profiles(): tipi senza distribuzione dei totali: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Trascritti per cellula (Step 2, passo 2)
#'
#' Ogni cellula riceve (1) un totale: quantile casuale della distribuzione dei trascritti per cellula del suo tipo (nuclei veri di R6,
#' D2-3) per il fattore di profondita' del tessuto (D2-9), arrotondato, almeno 1; (2) un sottotipo, estratto dalle quote del suo tipo
#' (D2-8); (3) i conteggi per gene: multinomiale (totale, profilo del sottotipo) con `noise = "M"`; con `noise = "G"` il profilo e' prima
#' moltiplicato gene per gene da un fattore gamma di media 1 e varianza `phi` e rinormalizzato (binomiale negativa per gene con totale
#' fissato). Le cellule si lavorano a blocchi di `chunk_size`; ogni blocco ha il suo flusso casuale L'Ecuyer-CMRG derivato da `seed`,
#' quindi il risultato non dipende da `workers`. Lo stato del generatore casuale dell'utente e' ripristinato all'uscita.
#' @param cell_df data.frame con `cell_id` e `cell_type` (uscita dello Step 1).
#' @param profiles uscita di `load_cell_profiles()`.
#' @param seed intero.
#' @param depth_factor fattore di profondita' (predefinito: quello dei profili, regola R-3 di S2.2).
#' @param noise `"M"` o `"G"` (predefinito: il vincitore della gara di S2.2 per il tessuto).
#' @param use_subtypes `FALSE` usa un profilo per tipo (media dei profili dei sottotipi pesata con le quote di cellule; la gara di S2.2
#'   pesava anche con i trascritti, quindi i due profili differiscono di poco); per le controprove.
#' @param totals_by `"type"`: totale dalla distribuzione del tipo; `"subtype"`: prima il sottotipo, poi il totale dalla distribuzione
#'   del sottotipo (`profiles$totals_subtype`; analisi PH-3 di S2.2, DEC-046). Predefinito (NULL): quello scritto nei metadati dei
#'   profili (`"subtype"` dove i profili vengono dai nuclei veri di R6, `"type"` per l'atlante). Con `use_subtypes = FALSE` vale sempre `"type"`.
#' @param chunk_size cellule per blocco (fa parte della definizione del risultato).
#' @param workers processi (`parallel::mclapply`).
#' @return lista `counts` (dgCMatrix geni x cellule, nomi `gene_id` e `cell_id`), `cells` (data.frame `cell_id`, `cell_type`, `subtype`,
#'   `n_transcripts`), `metadata`.
#' @export
simulate_cell_profiles <- function(cell_df, profiles, seed, depth_factor = profiles$depth_factor, noise = profiles$noise,
                                   use_subtypes = TRUE, totals_by = NULL, chunk_size = 2000L, workers = 1L) {
  if (is.null(totals_by)) totals_by <- if (!is.null(profiles$meta$totals_by) && !is.na(profiles$meta$totals_by)) profiles$meta$totals_by else "type"
  totals_by <- match.arg(totals_by, c("type", "subtype"))
  if (!is.data.frame(cell_df) || !all(c("cell_id", "cell_type") %in% names(cell_df)))
    stop("simulate_cell_profiles(): `cell_df` deve avere `cell_id` e `cell_type`.", call. = FALSE)
  if (anyNA(cell_df$cell_type)) stop("simulate_cell_profiles(): `cell_type` con NA.", call. = FALSE)
  if (anyNA(cell_df$cell_id) || anyDuplicated(cell_df$cell_id)) stop("simulate_cell_profiles(): `cell_id` con NA o duplicati.", call. = FALSE)
  if (missing(seed) || length(seed) != 1L || !is.finite(seed)) stop("simulate_cell_profiles(): serve un `seed`.", call. = FALSE)
  noise <- match.arg(noise, c("M", "G"))
  if (!is.numeric(depth_factor) || length(depth_factor) != 1L || !is.finite(depth_factor) || depth_factor <= 0)
    stop("simulate_cell_profiles(): `depth_factor` deve essere un numero > 0.", call. = FALSE)
  types <- unique(as.character(cell_df$cell_type))
  .s22_check_profiles(profiles, types)
  if (!use_subtypes) totals_by <- "type"
  if (totals_by == "subtype") {
    ts <- profiles$totals_subtype
    if (is.null(ts) || !identical(ts$q, profiles$totals$q) || length(setdiff(profiles$subtypes$subtype, names(ts))))
      stop("simulate_cell_profiles(): `totals_by = \"subtype\"` richiede i quantili dei totali per ogni sottotipo.", call. = FALSE)
  }
  P <- profiles$P; sub <- profiles$subtypes
  if (!use_subtypes) {                                         # un profilo per tipo: media dei sottotipi pesata con le quote
    tt <- unique(sub$cell_type)
    P <- t(vapply(tt, function(t) { k <- sub$cell_type == t; colSums(P[k, , drop = FALSE] * sub$share[k]) }, numeric(ncol(P))))
    P <- P / rowSums(P); rownames(P) <- tt
    sub <- data.frame(subtype = tt, cell_type = tt, share = 1, stringsAsFactors = FALSE)
  }
  G <- ncol(P); n <- nrow(cell_df)
  phi <- profiles$phi; gpos <- phi > 0; shape <- ifelse(gpos, 1 / ifelse(gpos, phi, 1), Inf)
  qs <- profiles$totals$q
  ctype <- as.character(cell_df$cell_type)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  # flussi casuali per blocco (L'Ecuyer-CMRG), indipendenti dal numero di processi
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  streams <- vector("list", length(chunks)); s <- get(".Random.seed", envir = globalenv())
  for (k in seq_along(chunks)) { streams[[k]] <- s; s <- parallel::nextRNGStream(s) }
  one_chunk <- function(k) {
    assign(".Random.seed", streams[[k]], envir = globalenv())
    idx <- chunks[[k]]; m <- length(idx); ty <- ctype[idx]
    # (1) totali: quantile casuale della distribuzione del tipo, per il fattore di profondita'; (2) sottotipi
    u <- stats::runif(m); Nc <- integer(m); si <- integer(m)
    draw_sub <- function() for (t in unique(ty)) { w <- which(ty == t); r <- which(sub$cell_type == t)
      si[w] <<- if (length(r) == 1L) r else r[sample.int(length(r), length(w), replace = TRUE, prob = sub$share[r])] }
    if (totals_by == "type") {
      for (t in unique(ty)) { w <- which(ty == t); Nc[w] <- as.integer(pmax(1, round(depth_factor * stats::approx(qs, profiles$totals[[t]], u[w])$y))) }
      draw_sub()
    } else {                                                     # prima il sottotipo, poi il totale dalla sua distribuzione
      draw_sub()
      for (s_ in unique(si)) { w <- which(si == s_); Nc[w] <- as.integer(pmax(1, round(depth_factor * stats::approx(qs, profiles$totals_subtype[[sub$subtype[s_]]], u[w])$y))) }
    }
    # (3) conteggi per gene
    ii <- vector("list", m); xx <- vector("list", m)
    if (noise == "M") {
      for (s_ in unique(si)) { w <- which(si == s_); Tt <- sum(Nc[w])
        draws <- sample.int(G, Tt, replace = TRUE, prob = P[s_, ])
        cell <- rep(seq_along(w), Nc[w])
        tab <- split(draws, cell)
        for (h in seq_along(w)) { tb <- tabulate(tab[[h]], nbins = G); nz <- which(tb > 0); ii[[w[h]]] <- nz; xx[[w[h]]] <- tb[nz] } }
    } else {
      for (h in seq_len(m)) { lam <- P[si[h], ]; z <- which(lam > 0 & gpos)
        lam[z] <- lam[z] * stats::rgamma(length(z), shape = shape[z], rate = shape[z]); lam <- lam / sum(lam)
        tb <- as.integer(stats::rmultinom(1L, Nc[h], lam)); nz <- which(tb > 0); ii[[h]] <- nz; xx[[h]] <- tb[nz] }
    }
    list(i = ii, x = xx, Nc = Nc, si = si)
  }
  out <- if (workers > 1L) parallel::mclapply(seq_along(chunks), one_chunk, mc.cores = workers, mc.preschedule = FALSE)
         else lapply(seq_along(chunks), one_chunk)
  bad <- vapply(out, function(o) inherits(o, "try-error"), logical(1))
  if (any(bad)) stop("simulate_cell_profiles(): errore in un blocco: ", as.character(out[[which(bad)[1]]]), call. = FALSE)
  ii <- unlist(lapply(out, `[[`, "i"), recursive = FALSE); xx <- unlist(lapply(out, `[[`, "x"), recursive = FALSE)
  Nc <- unlist(lapply(out, `[[`, "Nc")); si <- unlist(lapply(out, `[[`, "si"))
  lens <- lengths(ii)
  counts <- Matrix::sparseMatrix(i = unlist(ii), j = rep(seq_len(n), lens), x = as.numeric(unlist(xx)), dims = c(G, n),
                                 dimnames = list(profiles$genes$gene_id, as.character(cell_df$cell_id)))
  list(counts = counts,
       cells = data.frame(cell_id = cell_df$cell_id, cell_type = ctype, subtype = sub$subtype[si], n_transcripts = Nc, stringsAsFactors = FALSE),
       metadata = list(step = "S2.2", seed = seed, noise = noise, depth_factor = depth_factor, use_subtypes = use_subtypes, totals_by = totals_by,
                       chunk_size = chunk_size, archetype = profiles$archetype, profiles_winner = profiles$meta$winner))
}
