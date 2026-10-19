# Build thesis/data_inventory.csv from the CZ CELLxGENE Discover catalogue (week 1, T1).
#
# Usage (from thesis/):  Rscript scripts/00_data_inventory_cellxgene.R data_inventory.csv
# Downloads the public dataset catalogue (~10 MB) to a temp file; nothing is saved in the repo except the CSV.
# The selection (sel) and first-publication dates (first) were checked by hand on 2026-10-15; "TO VERIFY" = not yet confirmed.
# GEM-1 cutoff used: single-cell training = CELLxGENE as of July 2024 (GEM-1 preprint, doi:10.1101/2025.09.08.674753),
# to be confirmed by Synthesize Bio. Arm G held-out donors must be first public AFTER that date (PROTOCOLLO.md §9).
# Needs: jsonlite.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "data_inventory.csv"

cat_file <- file.path(tempdir(), "cellxgene_datasets.json")
if (!file.exists(cat_file)) {
  download.file("https://api.cellxgene.cziscience.com/curation/v1/datasets?visibility=PUBLIC",
                cat_file, mode = "wb", quiet = TRUE)
}
d <- jsonlite::fromJSON(cat_file, simplifyVector = FALSE)

rx <- c(colon    = "colon|rectum|caecum|cecum|large intestine|colorectal",
        lymphoid = "lymph node|tonsil",
        brain    = "cortex|gyrus|frontal|temporal|cerebral|brain")

# first CELLxGENE publication date, from /collections/{id}/versions (checked 2026-10-15)
first <- c("4d82bd8e" = "2025-08-05", "1cbfb478" = "2025-10-16", "3a844375" = "2026-06-08", "75e35754" = "2026-10-05",
           "64b24fda" = "2023-01-05", "33d19f34" = "2023-01-07", "7681c7d7" = "2022-06-15", "e5f58829" = "2021-07-14",
           "854c0855" = "2024-10-23", "69f13fde" = "2025-07-08", "131c68c1" = "2025-12-01", "cc431242" = "2025-06-06",
           "3a2af25b" = "2023-02-23", "283d65eb" = "2023-01-17", "d17249d2" = "2023-10-20", "35928d1c" = "2024-08-05",
           "f406a653" = "2025-01-31", "4624894f" = "2026-02-24", "7c7bd6c2" = "2025-11-10", "5c868b6f" = "2023-03-31")
V <- "TO VERIFY"

# one row per candidate: tissue, collection id prefix, name, original public date, GEM-1 status, role, notes
sel <- list(
  c("colon", "4d82bd8e", "Pelka 2021 (HTAN CRC)", "paper 2021 (Cell)", "PRE: public since 2021 - not for arm G held-out", "E/M training or E' (CRC + normal colon)", "Tumour + adjacent normal; on CELLxGENE only since 2025 but data public since 2021"),
  c("colon", "1cbfb478", "Moorman 2024 (CRC metastasis)", paste("paper online 2024-10-30 (Nature); preprint 2023-08; raw data deposit date", V), paste("POST? only if data first deposited after 2024-07 -", V), "candidate held-out for arm G (CRC + normal colon)", "Matched normal colon / primary / metastasis"),
  c("colon", "33d19f34", "Smillie 2019 (UC colon)", "paper 2019 (Cell)", "PRE", "E/M training (normal colon)", "Healthy donors of a UC study"),
  c("colon", "64b24fda", "Burclaff 2022 (healthy intestine)", "paper 2022 (CMGH)", "PRE", "E/M training or E' (normal colon)", "Proximal-to-distal survey"),
  c("colon", "7681c7d7", "James 2020 (colon immune)", "paper 2020 (Nat Immunol)", "PRE", "E/M training (normal colon)", "Immune-focused"),
  c("colon", "5c868b6f", "Kong 2023 (Crohn disease)", "paper 2023 (Immunity)", "PRE", "E/M training (normal colon)", "Use healthy donors only"),
  c("colon", "7c7bd6c2", "Garrido-Trigo 2023 (IBD + healthy)", "paper 2023 (Nat Commun)", "PRE", "E/M training (normal colon)", "Use healthy controls only"),
  c("colon", "75e35754", "Tearle 2026 (PSC/UC + healthy colon)", paste("paper 2026 (Nat Commun); medRxiv 2024-07-23; data deposit date", V), paste("POST? -", V), "candidate held-out for arm G (normal colon)", "Use healthy controls only"),
  c("colon", "3a844375", "Marteau 2026 (CRC atlas, integrated)", "integrates earlier public studies", "PRE (re-uses older data)", "avoid as held-out", "Integrated atlas: donors overlap other studies (e.g. Pelka) - leakage risk"),
  c("lymphoid", "e5f58829", "Tabula Sapiens v1 (lymph node)", "paper 2022 (Science); CELLxGENE 2021", "PRE: Tabula Sapiens v1 is in GEM-1 training - NEVER for arm G", "E/M training (lymph node)", "Named as GEM-1 single-cell training data in the preprint"),
  c("lymphoid", "854c0855", "Xu 2023 CellHint (lymph node, integrated)", "integrates earlier public studies", "PRE (re-uses older data)", "avoid as held-out", "Integrated: overlaps Tabula Sapiens and others"),
  c("lymphoid", "3a2af25b", "King 2021 (tonsil)", "paper 2021 (Sci Immunol)", "PRE", "E/M training (tonsil)", "Tonsil B-cell maturation"),
  c("lymphoid", "69f13fde", "Fergusson 2025 LNCS (healthy lymph nodes)", paste("paper 2025 (Cell Rep); data deposit date", V), paste("POST? -", V), "candidate held-out for arm G (lymph node)", "Slices cultured with/without adjuvant: use unstimulated only"),
  c("lymphoid", "131c68c1", "Siu 2025 (axillary lymph node FNA)", paste("paper 2025 (eBioMedicine); data deposit date", V), paste("POST? -", V), "candidate held-out for arm G (lymph node)", "Fine-needle aspirates: few stromal cells"),
  c("lymphoid", "cc431242", "Wells 2025 (immune aging, multi-tissue)", paste("paper 2025 (Nat Immunol); data deposit date", V), paste("POST? -", V), "candidate (lymph node subset)", "Check which donors include lymph node"),
  c("brain", "283d65eb", "Siletti 2023 (Human Brain Cell Atlas)", "paper 2023 (Science)", "PRE", "E/M training (adult cortex)", "Few donors, many regions"),
  c("brain", "d17249d2", "Jorstad 2023 (cortical areas)", "paper 2023 (Science)", "PRE", "E/M training or E' (adult cortex)", ""),
  c("brain", "35928d1c", "Johansen 2023 (Human variation study)", "paper 2023 (Science); CELLxGENE only 2024-08-05", "PRE: data public with the 2023 paper", "E/M training or E' (adult neocortex)", "Neurosurgical neocortex, many donors"),
  c("brain", "f406a653", "Clarence 2025 (postnatal brain development)", paste("paper 2025 (Nat Genet); data deposit date", V), paste("POST? -", V), "candidate held-out for arm G (cortex)", "Postnatal ages: keep adult donors only?"),
  c("brain", "4624894f", "Galvao 2026 (adolescent brain, bioRxiv)", paste("preprint 2026-02; deposit date", V), "POST (preprint 2026) - confirm deposit", "candidate held-out for arm G (cortex)", "Normal tissue from surgery; adolescents; not peer reviewed")
)

labels <- function(x) vapply(x, function(o) o$label, character(1))
csort  <- function(x) sort(unique(x), method = "radix")   # locale-independent order
coll   <- substr(vapply(d, function(e) e$collection_id, character(1)), 1, 8)

rows <- lapply(sel, function(s) {
  g <- s[1]; p <- s[2]
  # datasets of this collection: human, with a matching tissue (colon also: matching disease)
  es <- Filter(function(e) {
    "Homo sapiens" %in% labels(e$organism) &&
      (any(grepl(rx[[g]], labels(e$tissue), ignore.case = TRUE)) ||
         (g == "colon" && any(grepl("colorectal|colon", labels(e$disease), ignore.case = TRUE))))
  }, d[coll == p])
  don   <- unique(unlist(lapply(es, function(e) unlist(e$donor_id))))
  ct    <- unique(unlist(lapply(es, function(e) labels(e$cell_type))))
  ass   <- csort(unlist(lapply(es, function(e) labels(e$assay))))
  dis   <- csort(unlist(lapply(es, function(e) labels(e$disease))))
  cells <- sum(vapply(es, function(e) as.numeric(e$cell_count), numeric(1)))
  e0 <- es[[1]]
  c(tissue_group = g, dataset = s[3],
    publication = paste0("https://doi.org/", if (is.null(e0$collection_doi)) "" else e0$collection_doi),
    source_url = paste0("https://cellxgene.cziscience.com/collections/", e0$collection_id),
    donors_n = length(don), cells_n_max = format(cells, scientific = FALSE), cell_types_n_max = length(ct),
    technology = paste(ass, collapse = "; "), disease_in_data = paste(dis, collapse = "; "),
    first_public_cellxgene = first[[p]], original_public_date_evidence = s[4], gem1_cutoff_status = s[5],
    licence = "CC BY 4.0 (CELLxGENE)", candidate_role = s[6],
    notes = paste0(if (nzchar(s[7])) paste0(s[7], " | ") else "",
                   "counts summed over matching CELLxGENE datasets (subsets may overlap); all donors, any disease"),
    checked_on = "2026-10-15")
})
cols <- names(rows[[1]])

# two candidates not taken from the CELLxGENE catalogue
extra <- list(
  c("colon", "Oliveira 2025 / 10x (CRC Flex, matched to Visium HD CRC)",
    "https://doi.org/10.1101/2024.06.04.597233 (preprint; Nat Genet 2025)", "GEO GSE280311 (part of GSE280318)",
    V, V, V, "10x Chromium Flex (FFPE, probe-based)", "CRC; normal adjacent", "not on CELLxGENE",
    paste("GEO release 2025-03-18; 2024-06 preprint pointed to data on 10x website - first public date", V),
    paste("UNCLEAR: GEO after cutoff, 10x website possibly before -", V), V,
    "candidate; same patients as the main project's Visium HD CRC", "Flex is probe-based (fixed gene panel)", "2026-10-15"),
  c("lymphoid", "Massoni-Badosa 2024 (tonsil atlas, HCATonsilData)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10869140",
    "https://bioconductor.org/packages/HCATonsilData", "17 (6 children, 8 young adults, 3 older adults)",
    paste(V, "(paper versions differ)"), V, "10x scRNA-seq (+ multiome, CITE-seq, scATAC)", "normal tonsil",
    "not found on CELLxGENE", "paper 2024-02-13 (Immunity); preprint 2022-06", "PRE", V,
    "E/M training (tonsil)", "Rich annotation with a cell-type glossary", "2026-10-15")
)
extra <- lapply(extra, setNames, cols)

inv <- as.data.frame(do.call(rbind, c(rows, extra)), stringsAsFactors = FALSE)
write.csv(inv, out, row.names = FALSE, fileEncoding = "UTF-8")
cat("ok:", nrow(inv), "rows ->", out, "\n")
