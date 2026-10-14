# tools/S1.3d_f4.R — S1.3d (BL-084), candidato F4: taglio del Voronoi in un colpo solo (cs = 0, una regione).
# Pre-registrazione: results/S1.3d/S1.3d_preregistration.md (ee6208b). Stessi tile di R/04b3 (.tv_window, .tv_tiles);
# lati di tutti i tile + bordo della regione annodati UNA volta (st_union delle linee) e poligonizzati: ogni punto di
# intersezione lato/bordo è calcolato una volta sola, quindi i lati condivisi hanno vertici identici per costruzione.
# Frammenti come D-S1.3.2 / .tv_fragments: il pezzo con il generatore resta; gli orfani vanno, a passate simultanee, alla
# cella con il confine condiviso più lungo (> 1e-6 µm; parità -> indice minore); senza confinante restano al generatore.
# Confine condiviso = somma delle lunghezze dei lati comuni fra facce (lati identici: stessa rete).
f4_build <- function(x, y, reg, lmin = 1e-6) {
  t0 <- proc.time()[[3]]; n <- length(x)
  rw <- .tv_window(reg); tl <- .tv_tiles(x, y, rw)
  lines <- c(st_cast(st_boundary(tl), "MULTILINESTRING"), st_cast(st_boundary(reg), "MULTILINESTRING"))
  net <- st_union(lines)                                            # annodamento unico (OverlayNG, virgola mobile)
  faces <- st_collection_extract(st_polygonize(net), "POLYGON")
  faces <- st_cast(faces, "POLYGON")
  fp <- st_point_on_surface(faces)
  in_reg <- lengths(st_intersects(fp, reg)) > 0
  faces <- faces[in_reg]; fp <- fp[in_reg]; nf <- length(faces)
  hit <- st_intersects(fp, tl); nh <- lengths(hit)
  owner <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, 1L)
  amb <- which(nh != 1)
  for (f in amb) {                                                  # punto interno ambiguo: tile con la massima area comune
    cand <- unlist(st_intersects(faces[f], tl))
    owner[f] <- cand[which.max(vapply(cand, function(k) sum(as.numeric(st_area(st_intersection(faces[f], tl[k])))), 0))]
  }
  fa <- as.numeric(st_area(faces))
  # adiacenza fra facce lungo lati (interni disgiunti, bordi con intersezione 1-dimensionale) e lunghezza condivisa
  adj <- st_relate(faces, faces, pattern = "F***1****")
  P <- do.call(rbind, lapply(seq_along(adj), function(a) { b <- adj[[a]]; b <- b[b > a]; if (length(b)) cbind(a, b) }))
  if (is.null(P)) P <- matrix(integer(0), 0, 2)
  fb <- st_boundary(faces)
  L <- if (nrow(P)) vapply(seq_len(nrow(P)), function(q) { s <- st_intersection(fb[P[q, 1]], fb[P[q, 2]]); if (!length(s)) 0 else sum(as.numeric(st_length(s))) }, 0) else numeric(0)
  # componenti connesse delle facce di ogni cella
  par <- seq_len(nf); fnd <- function(a) { while (par[a] != a) { par[a] <<- par[par[a]]; a <- par[a] }; a }
  same <- if (nrow(P)) owner[P[, 1]] == owner[P[, 2]] else logical(0)
  for (q in which(same)) { ra <- fnd(P[q, 1]); rb <- fnd(P[q, 2]); if (ra != rb) par[max(ra, rb)] <- min(ra, rb) }
  comp <- vapply(seq_len(nf), fnd, 1L)
  cell_of <- owner                                                 # cella attuale di ogni faccia
  gen <- st_sfc(lapply(seq_len(n), function(i) st_point(c(x[i], y[i]))))
  orphans <- list(); lost <- integer(n); gained <- integer(n)
  for (i in unique(owner)) {
    fi <- which(owner == i); cs <- unique(comp[fi]); if (length(cs) < 2) next
    gin <- unlist(st_intersects(gen[i], faces[fi])); own <- if (length(gin)) comp[fi[gin[1]]] else {
      ca <- vapply(cs, function(cc) sum(fa[fi[comp[fi] == cc]]), 0); cs[which.max(ca)] }
    for (cc in setdiff(cs, own)) orphans[[length(orphans) + 1L]] <- list(f = fi[comp[fi] == cc], from = i)
    lost[i] <- length(cs) - 1L
  }
  for (o in seq_along(orphans)) cell_of[orphans[[o]]$f] <- NA_integer_   # come .tv_fragments: orfani in attesa senza cella
  nfrag <- length(orphans); n_re <- 0L; margins <- numeric(0); isolated <- integer(0); pending <- seq_len(nfrag)
  if (nfrag) for (pass in 1:50) {
    dec <- vapply(pending, function(o) {
      ff <- orphans[[o]]$f; from <- orphans[[o]]$from
      q <- which((P[, 1] %in% ff) != (P[, 2] %in% ff))             # lati fra l'orfano e il resto
      if (!length(q)) return(NA_integer_)
      other <- ifelse(P[q, 1] %in% ff, P[q, 2], P[q, 1]); oc <- cell_of[other]
      len <- tapply(L[q], oc, sum); ids <- as.integer(names(len)); keep <- ids != from
      len <- len[keep]; ids <- ids[keep]; ok <- len > lmin
      if (!any(ok)) return(NA_integer_)
      best <- ids[ok & len == max(len[ok])][1]
      srt <- sort(len[ok], decreasing = TRUE); margins[[length(margins) + 1L]] <<- if (length(srt) > 1) srt[1] - srt[2] else Inf
      best }, 1L)
    if (all(is.na(dec))) break
    for (k in which(!is.na(dec))) { o <- pending[k]; cell_of[orphans[[o]]$f] <- dec[k]; gained[dec[k]] <- gained[dec[k]] + 1L; n_re <- n_re + 1L }
    pending <- pending[is.na(dec)]; if (!length(pending)) break
  }
  for (o in pending) { cell_of[orphans[[o]]$f] <- orphans[[o]]$from; lost[orphans[[o]]$from] <- lost[orphans[[o]]$from] - 1L; isolated <- c(isolated, orphans[[o]]$from) }
  geoms <- lapply(seq_len(n), function(i) {
    fi <- which(cell_of == i); if (!length(fi)) return(st_polygon())
    if (length(fi) == 1) return(faces[[fi]])
    u <- st_union(faces[fi], is_coverage = TRUE)[[1]]
    if (inherits(u, "MULTIPOLYGON") && length(unclass(u)) == 1) st_polygon(unclass(u)[[1]]) else u
  })
  out <- st_sfc(geoms)
  list(terr = out, info = list(n_faces = nf, n_ambiguous = length(amb), n_fragments = nfrag, n_reassigned = n_re,
       n_isolated = length(unique(isolated)), min_margin = if (length(margins)) min(margins) else Inf, n_tie = sum(margins < lmin),
       n_invalid = sum(!st_is_valid(out)), n_empty = sum(st_is_empty(out)), lost = lost, gained = gained, elapsed_s = proc.time()[[3]] - t0))
}
