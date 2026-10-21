# tools/S1.5_faults.R - S1.5: guasti costruiti per i controlli del livello cellulare (test G-1 e banco G-4).
# Ogni guasto modifica un cell_layer valido in un punto e con un'area nota delta (µm²); restituisce il cell_layer guasto,
# il controllo atteso e la cellula toccata. Quadrati di lato sqrt(delta) posti con st_point_on_surface().

.f_sq <- function(p, delta) { h <- sqrt(delta) / 2; xy <- sf::st_coordinates(p)[1, 1:2]
  sf::st_polygon(list(rbind(xy + c(-h, -h), xy + c(h, -h), xy + c(h, h), xy + c(-h, h), xy + c(-h, -h)))) }
.f_set <- function(sfc, i, g) { sfc[[i]] <- g; sf::st_sfc(sfc) }
.f_union <- function(a, b) sf::st_union(sf::st_sfc(a), sf::st_sfc(b))[[1]]
.f_diff <- function(a, b) sf::st_difference(sf::st_sfc(a), sf::st_sfc(b))[[1]]

# coppia (i, j) di territori interni vicini, con j abbastanza grande per contenere il quadrato
.f_pair <- function(cl, delta) {
  td <- cl$territory_df; ok <- which(!td$clipped & td$n_pieces_gained == 0 & td$territory_area > 20 * delta + 1)
  nb <- sf::st_touches(cl$cell_territories[ok], cl$cell_territories)
  for (k in seq_along(ok)) { j <- intersect(nb[[k]], ok); if (length(j)) return(c(i = ok[k], j = j[1])) }
  stop("nessuna coppia adatta")
}

#' Costruisce un guasto di tipo `type` e area `delta` su un cell_layer valido
make_fault <- function(cl, type = c("overlap", "gap", "outside_region", "nucleus_out", "nuclei_overlap", "nan_area", "invalid"),
                       delta = 0.1) {
  type <- match.arg(type); Tt <- cl$cell_territories; Nn <- cl$cell_nuclei
  if (type %in% c("overlap", "gap", "nucleus_out", "nuclei_overlap", "nan_area")) pij <- .f_pair(cl, delta)
  if (type == "overlap") {        # T_i invade T_j per delta
    j <- pij[["j"]]; S <- .f_sq(sf::st_point_on_surface(sf::st_sfc(.f_diff(Tt[[j]], Nn[[j]]))), delta)
    cl$cell_territories <- .f_set(Tt, pij[["i"]], .f_union(Tt[[pij[["i"]]]], S)); return(list(cl = cl, expect = "C2", cell = pij[["i"]]))
  }
  if (type == "gap") {            # buco di area delta dentro T_j, fuori dal nucleo
    j <- pij[["j"]]; S <- .f_sq(sf::st_point_on_surface(sf::st_sfc(.f_diff(Tt[[j]], sf::st_buffer(Nn[j], 1)[[1]]))), delta)
    cl$cell_territories <- .f_set(Tt, j, .f_diff(Tt[[j]], S)); return(list(cl = cl, expect = "C2", cell = j))
  }
  if (type == "outside_region") { # un territorio di bordo esce dalla regione per delta
    td <- cl$territory_df; i <- which(td$clipped)[1]
    R <- cl$region_polygons[match(cl$cell_df$region_id[i], cl$region_df$region_id)]
    out <- sf::st_difference(sf::st_buffer(Tt[i], 2), R)
    S <- .f_sq(sf::st_point_on_surface(sf::st_sfc(sf::st_difference(out, sf::st_buffer(R, 0.5))[[1]])), delta)
    cl$cell_territories <- .f_set(Tt, i, .f_union(Tt[[i]], S)); return(list(cl = cl, expect = "C2", cell = i))
  }
  if (type == "nucleus_out") {    # il nucleo di i esce nel territorio di j, fuori dal nucleo di j
    j <- pij[["j"]]; S <- .f_sq(sf::st_point_on_surface(sf::st_sfc(.f_diff(Tt[[j]], sf::st_buffer(Nn[j], 1)[[1]]))), delta)
    cl$cell_nuclei <- .f_set(Nn, pij[["i"]], .f_union(Nn[[pij[["i"]]]], S)); return(list(cl = cl, expect = "C3", cell = pij[["i"]]))
  }
  if (type == "nuclei_overlap") { # il nucleo di i invade il nucleo di j
    j <- pij[["j"]]; S <- .f_sq(sf::st_point_on_surface(Nn[j]), delta)
    cl$cell_nuclei <- .f_set(Nn, pij[["i"]], .f_union(Nn[[pij[["i"]]]], S)); return(list(cl = cl, expect = "C3", cell = pij[["i"]]))
  }
  if (type == "nan_area") { cl$cell_df$nucleus_area[pij[["i"]]] <- NaN; return(list(cl = cl, expect = "C5", cell = pij[["i"]])) }
  # invalid: il territorio di i diventa un poligono a farfalla (auto-intersecante) sul suo riquadro
  i <- which(!cl$territory_df$clipped)[1]; b <- sf::st_bbox(Tt[i])
  bow <- sf::st_polygon(list(rbind(c(b$xmin, b$ymin), c(b$xmax, b$ymax), c(b$xmax, b$ymin), c(b$xmin, b$ymax), c(b$xmin, b$ymin))))
  cl$cell_territories <- .f_set(Tt, i, bow); list(cl = cl, expect = "C9", cell = i)
}

FAULT_TYPES <- c("overlap", "gap", "outside_region", "nucleus_out", "nuclei_overlap", "nan_area", "invalid")
