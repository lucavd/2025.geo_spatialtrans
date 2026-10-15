# tools/S1.3d_compare.R — S1.3d (BL-084): confronto delle uscite dei candidati con la verità esatta (tools/S1.3d_arbiter.py).
# Tolleranze pre-registrate (g = 1e-6 per tutti): unione/diff. simmetrica tau_U = g*(P_reg + sum P); coppia tau_ij = g*(P_i + P_j);
# cella tau_i = g*P_i; pezzo tau_kb = g*perimetro del pezzo esatto (anelli S-H); pezzo spurio: g*perimetro del candidato (F6: g*4s);
# O3 tau = somma dei tau dei sotto-pezzi esatti a 2 µm (approssimata con il tau del pezzo esatto da 8 µm + perimetri interni: vedi dev_notes).
#' O1: restituisce list(summary = data.frame di 1 riga, viol = data.frame delle violazioni)
cmp_o1 <- function(o, tr, tp, tc, P, Preg, g = G_GRID) {
  tauU <- g * (Preg + sum(P)); tij <- g * (P[tp$i] + P[tp$j]); ti <- g * P
  eU <- abs(o$union_a - tr$union); eS <- abs(o$symd_a - tr$symd)
  eij <- abs(o$pij - tp$a); eji <- abs(o$pji - tp$a); ec <- abs(o$area - tc$a)
  vij <- which(is.na(eij) | eij > tij); vji <- which(is.na(eji) | eji > tij); vc <- which(is.na(ec) | ec > ti)
  exc <- length(o$errors) + o$n_pair_err
  s <- data.frame(cand = o$cand, union_err = eU, union_viol = is.na(eU) || eU > tauU, symd_err = eS, symd_viol = is.na(eS) || eS > tauU, tauU = tauU,
                  n_pairs = nrow(tp), pair_viol_ij = length(vij), pair_viol_ji = length(vji), pair_viol_any = length(union(vij, vji)),
                  pair_max_ratio = max(c(eij / tij, eji / tij), na.rm = TRUE), cell_viol = length(vc), cell_max_ratio = max(ec / ti, na.rm = TRUE),
                  n_exceptions = exc, t_prepare = o$time[["prepare"]], t_union = o$time[["union"]], t_pairs = o$time[["pairs"]], stringsAsFactors = FALSE)
  s$silent <- (s$union_viol || s$symd_viol || s$pair_viol_any > 0 || s$cell_viol > 0) && exc == 0
  s$fail <- s$union_viol || s$symd_viol || s$pair_viol_any > 0 || s$cell_viol > 0 || exc > 0
  viol <- rbind(if (length(vij)) data.frame(what = "pair_ij", i = tp$i[vij], j = tp$j[vij], cand_val = o$pij[vij], truth = tp$a[vij], tau = tij[vij]),
                if (length(vji)) data.frame(what = "pair_ji", i = tp$j[vji], j = tp$i[vji], cand_val = o$pji[vji], truth = tp$a[vji], tau = tij[vji]),
                if (length(vc)) data.frame(what = "cell", i = vc, j = NA, cand_val = o$area[vc], truth = tc$a[vc], tau = ti[vc]),
                if (isTRUE(s$union_viol)) data.frame(what = "union", i = NA, j = NA, cand_val = o$union_a, truth = tr$union, tau = tauU),
                if (isTRUE(s$symd_viol)) data.frame(what = "symd", i = NA, j = NA, cand_val = o$symd_a, truth = tr$symd, tau = tauU))
  list(summary = s, viol = viol)
}
#' O2/O3: pezzi del candidato contro pezzi esatti (chiave i, ix, iy)
cmp_o2 <- function(pc, tr, s, cand, g = G_GRID, tau_truth = NULL) {
  if (is.null(pc)) return(list(summary = data.frame(cand = cand, s = s, n_truth = nrow(tr), n_cand = NA, n_viol = NA, n_missing = NA, n_spurious = NA,
                                                     max_ratio = NA, sum_signed = NA, exception = TRUE, fail = TRUE, silent = FALSE), viol = NULL))
  kt <- sprintf("%d %d %d", as.integer(tr$i), as.integer(tr$ix), as.integer(tr$iy)); kc <- sprintf("%d %d %d", as.integer(pc$i), as.integer(pc$ix), as.integer(pc$iy))   # interi: paste(100000) = "1e+05" (bug corretto, D-3)
  if (anyDuplicated(kc)) { ag <- tapply(pc$area, kc, sum); pp <- tapply(pc$perim, kc, function(v) if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE))
    pc <- data.frame(k = names(ag), area = as.numeric(ag), perim = as.numeric(pp[names(ag)])) } else pc <- data.frame(k = kc, area = pc$area, perim = pc$perim)
  tt <- data.frame(k = kt, area_t = tr$area, perim_t = tr$perim)
  m <- merge(tt, pc, by = "k", all = TRUE)
  m$area[is.na(m$area)] <- 0; ta <- m$area_t; ta[is.na(ta)] <- 0
  tau <- if (!is.null(tau_truth)) tau_truth[match(m$k, kt)] else g * m$perim_t
  spur <- is.na(m$area_t); tau[spur] <- g * ifelse(is.na(m$perim[spur]), 4 * s, m$perim[spur])
  err <- abs(m$area - ta); v <- err > tau
  sm <- data.frame(cand = cand, s = s, n_truth = nrow(tr), n_cand = nrow(pc), n_viol = sum(v), n_missing = sum(v & m$area == 0 & !spur),
                   n_spurious = sum(v & spur), max_ratio = max(err / tau), sum_signed = sum(m$area) - sum(ta), exception = FALSE, stringsAsFactors = FALSE)
  sm$fail <- sm$n_viol > 0; sm$silent <- sm$fail
  list(summary = sm, viol = if (any(v)) head(m[v, ], 200) else NULL)
}
