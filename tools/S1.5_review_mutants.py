#!/usr/bin/env python3
# tools/S1.5_review_mutants.py - S1.5, controprova dopo la revisione: i 19 mutanti del revisore del codice
# (review_code/make_mutants.py e make_r19.py), stesse modifiche sul codice corretto; R15 adattato alla nuova riga di C6
# (conteggio dalle cellule, revisione FND-04). Alberi in /mnt/micron/geo_spatialtrans/S1.5/postreview/tree.

import os, shutil, sys, json
REPO = "/home/user/2025.geo_spatialtrans"; RC = "/mnt/micron/geo_spatialtrans/S1.5/postreview"
V = "R/04b5_cell_layer_validation.R"; P = "R/04b_cell_placement.R"; T = "R/04b0_tissue_presets.R"; E = "R/04b6_estimate_pixel_size.R"
MUT = {
 "R01": (V, "nb_c3 <- .nbad(dn <= g * pN + g^2) + .nbad(wi)", "nb_c3 <- .nbad(dn <= g * pN + g^2)", "C3: st_within indipendente tolto"),
 "R02": (V, "nb_out <- .nbad(ao <= g * pT[bnd] + g^2)", "nb_out <- .nbad(ao <= 1e6 * (g * pT[bnd] + g^2))", "C2: tolleranza dell'uscita dalla regione x1e6"),
 "R03": (V, "s <- sf::st_intersects(g)", "s <- sf::st_touches(g)", "C2/C3: coppie con st_touches (esclude le coppie sovrapposte)"),
 "R04": (V, "nb_sum <- .nbad(rel_sum <= sum_rel_tol)", "nb_sum <- 0L", "C2: controllo della somma per regione spento"),
 "R05": (V, ".nbad <- function(ok) sum(is.na(ok) | !ok)", ".nbad <- function(ok) sum(!ok, na.rm = TRUE)", "NA/NaN contati come superati"),
 "R06": (P, 'u <- .pc_with_seed(ss[["pale"]], stats::runif(nrow(cd)))', 'u <- .pc_with_seed(ss[["seed"]], stats::runif(nrow(cd)))', "nuclei pallidi estratti con il seme della semina"),
 "R07": (T, "nucleus_area_slope = p$beta, nucleus_area_sdlog = p$sigma", "nucleus_area_slope = p$sigma, nucleus_area_sdlog = p$beta", "catalogo: beta e sigma scambiati (tabella fonti intatta)"),
 "R08": (T, "min_dist_um = p$min_dist_um,", "min_dist_um = 2 * p$min_dist_um,", "catalogo: distanza minima raddoppiata (tabella fonti intatta)"),
 "R09": (P, "cd$is_minority_type <- as.character(cd$cell_type) != rdf$dominant_cell_type[match(cd$region_id, rdf$region_id)]", "cd$is_minority_type <- rep(FALSE, nrow(cd))", "is_minority_type sempre FALSE"),
 "R10": (P, "tb <- tb[order(tb$comp, -tb$n, tb$code), ]", "tb <- tb[order(tb$comp, tb$n, tb$code), ]", "merge: cluster con MENO lati condivisi"),
 "R11": (P, 'if (!is.null(attr(clust, "stride"))) stride <- as.integer(attr(clust, "stride")) else {', 'if (FALSE) stride <- 1L else {', "BL-042: attributo stride ignorato"),
 "R12": (E, "if (is.finite(chk) && abs(chk / mpp - 1) > 1e-3)", "if (FALSE)", "metadati: controllo sulla griglia dei bin tolto"),
 "R13": (V, "nb <- sum(rel > 0.2 | !is.finite(rel), na.rm = TRUE)", "nb <- sum(rel > 2 | !is.finite(rel), na.rm = TRUE)", "C1: soglia 0.20 -> 2"),
 "R14": (P, 'size_conversion = c("preserve_area", "none"), validate', 'size_conversion = c("none", "preserve_area"), validate', "place_cells: braccio (b) come predefinito"),
 "R15": (V, "nb6 <- sum(rd$area_um2 >= min_region_area_um2 & ncr6 == 0)", "nb6 <- 0L", "C6 spento"),
 "R16": (E, 'if (N < 30 && .ep_mut() != "M9")', 'if (N < 3 && .ep_mut() != "M9")', "stima: soglia di rifiuto N < 30 -> N < 3"),
 "R17": (T, "p$rho * (p$density / p$density_model)^((1 - p$beta) / 2)", "p$rho * (p$density_model / p$density)^((1 - p$beta) / 2)", "rho': rapporto delle densita' invertito"),
 "R18": (V, "nb_ov <- 0L; max_ov <- 0\n  if (nrow(pn)) {", "nb_ov <- 0L; max_ov <- 0\n  if (FALSE) {", "C3: disgiunzione dei nuclei spenta"),
}
def make(tag, edit=None):
    d = os.path.join(RC, "tree", tag)
    if os.path.exists(d): shutil.rmtree(d)
    os.makedirs(d)
    shutil.copytree(os.path.join(REPO, "R"), os.path.join(d, "R"))
    os.makedirs(os.path.join(d, "tools")); shutil.copy(os.path.join(REPO, "tools/S1.5_faults.R"), os.path.join(d, "tools"))
    os.makedirs(os.path.join(d, "docs")); shutil.copy(os.path.join(REPO, "docs/BIO_REFERENCES.md"), os.path.join(d, "docs"))
    for sub, f in (("results/S1.4", "S1.4_calibration.csv"), ("results/R3b", "R3b_refs.csv")):
        os.makedirs(os.path.join(d, sub), exist_ok=True); shutil.copy(os.path.join(REPO, sub, f), os.path.join(d, sub))
    os.symlink(os.path.join(REPO, "renv"), os.path.join(d, "renv"))
    if edit:
        f, old, new, _ = edit; p = os.path.join(d, f); s = open(p).read()
        n = s.count(old); assert n == 1, (tag, n)
        open(p, "w").write(s.replace(old, new))
    return d
MUT["R19"] = ("R/04b3_tessellate_voronoi.R", "  terr <- fin$sfc\n", "  terr <- fin$sfc; terr[[1]] <- sf::st_buffer(terr[1], -0.2)[[1]]\n", "tassellazione: il territorio 1 perde una corona di 0.2 um (lacuna nella partizione)")
os.makedirs(RC, exist_ok=True)
make("base")
for k, e in MUT.items(): make(k, e)
json.dump({k: v[3] for k, v in MUT.items()}, open(os.path.join(RC, "mutants.json"), "w"), ensure_ascii=False, indent=1)
print("ok", len(MUT))
