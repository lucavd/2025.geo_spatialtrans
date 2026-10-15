# tools/S1.4c_decision.R — S1.4c: funzioni di verdetto, punteggio, bootstrap sui semi, AUC e CP-2 (pure, testate in
# R/testing/test_S1.4c_decision.R su dati sintetici prima del run). Regole dalla pre-registrazione (commit c41cf8b).
S14C_PANEL <- c("B1", "B2", "B7", "B8", "S1", "S2", "S3")

#' Verdetto di fold: distanza (mediana sui semi) contro il metro inter-ROI del ROI di test.
fold_verdict <- function(d, yard) ifelse(d <= stats::median(yard), "PASS", ifelse(d <= max(yard), "WARN", "FAIL"))
#' Verdetto d'archetipo dai 5 verdetti di fold (regola di S1.4): PASS se >= 4 PASS, FAIL se >= 2 FAIL, altrimenti WARN.
arch_verdict <- function(v) if (sum(v == "PASS") >= 4) "PASS" else if (sum(v == "FAIL") >= 2) "FAIL" else "WARN"

#' Verdetti per archetipo x metrica x variante.
#' @param dist data.frame(archetype, fold, seed_idx, variant, <metriche>) ; @param yards list[[A]][[fold]][[metrica]] = vettore
#' @param seeds_of opzionale: function(A, fold) -> indici dei semi da usare (bootstrap); default tutti.
verdicts <- function(dist, yards, metrics = S14C_PANEL, seeds_of = NULL) {
  out <- list()
  for (A in sort(unique(dist$archetype))) for (V in sort(unique(dist$variant))) for (m in metrics) {
    fv <- vapply(1:5, function(k) {
      z <- dist[dist$archetype == A & dist$variant == V & dist$fold == k, c("seed_idx", m)]
      z <- z[order(z$seed_idx), ]
      sel <- if (is.null(seeds_of)) z[[m]] else z[[m]][match(seeds_of(A, k), z$seed_idx)]
      fold_verdict(stats::median(sel), yards[[A]][[k]][[m]])
    }, "")
    out[[length(out) + 1]] <- data.frame(archetype = A, variant = V, metric = m, verdict = arch_verdict(fv),
                                         n_pass = sum(fv == "PASS"), n_warn = sum(fv == "WARN"), n_fail = sum(fv == "FAIL"))
  }
  do.call(rbind, out)
}
score_of <- function(vt, metrics) tapply(vt$verdict[vt$metric %in% metrics] == "PASS", vt$variant[vt$metric %in% metrics], sum)

#' Bootstrap sui semi: per replica, per (archetipo, fold) 10 semi con reinserimento, gli stessi per tutte le varianti
#' (le varianti condividono il flusso z del seme: confronto appaiato). @return matrice repliche x varianti dei punteggi.
boot_scores <- function(dist, yards, metrics, B = 200L, seed = 20261015L, nseed = 10L) {
  set.seed(seed); vs <- sort(unique(dist$variant))
  m <- t(vapply(seq_len(B), function(b) {
    draw <- list(); for (A in sort(unique(dist$archetype))) for (k in 1:5) draw[[paste(A, k)]] <- sample.int(nseed, nseed, replace = TRUE)
    sc <- score_of(verdicts(dist, yards, metrics, seeds_of = function(A, k) draw[[paste(A, k)]]), metrics)
    as.numeric(sc[vs])
  }, numeric(length(vs))))
  colnames(m) <- vs; m
}

#' Regola di decisione (punti 3-5): restituisce la variante scelta e la tabella dei criteri.
decide <- function(score, pboot, eligible, margin = 3, pmin = 0.9, order_simple = c("K-c", "K-1", "K-2")) {
  cand <- setdiff(order_simple, "K-c")
  tab <- data.frame(variant = cand, eligible = eligible[cand], score = as.numeric(score[cand]), score_Kc = as.numeric(score[["K-c"]]),
                    margin_ok = as.numeric(score[cand]) >= as.numeric(score[["K-c"]]) + margin,
                    p_boot = as.numeric(pboot[cand]), p_ok = as.numeric(pboot[cand]) >= pmin)
  tab$wins <- tab$eligible & tab$margin_ok & tab$p_ok
  w <- tab[tab$wins, ]
  choice <- if (!nrow(w)) "K-c" else { w <- w[order(-w$score, match(w$variant, order_simple)), ]; w$variant[1] }
  list(choice = choice, table = tab)
}

#' AUC (Mann-Whitney) dello score per l'esito binario.
auc <- function(score, y) { y <- as.logical(y); r <- rank(score); n1 <- sum(y); n0 <- sum(!y)
  if (!n1 || !n0) return(NA_real_); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

#' CP-2: per archetipo, numero di ROI in cui la distanza fra segmentatori supera il massimo del metro inter-ROI.
#' @param dseg vettore (5 ROI) ; @param yards list (5) di vettori D(primario k, primario j). Supera se >= 3/5 ROI.
cp2_exceeds <- function(dseg, yards) sum(dseg > vapply(yards, max, 0)) >= 3
