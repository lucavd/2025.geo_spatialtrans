# geo.spatialtrans: a Visium HD simulator with a cell layer

`geo.spatialtrans` simulates Visium HD spatial transcriptomics data in order to test analysis methods against a known truth.
Its distinctive feature is an explicit **cell layer** between the tissue image and the 2 µm capture grid: every simulated cell has a
position, a territory (a Voronoi tile clipped to its tissue region) and a nucleus. Step 2 generates the transcripts of each cell, moves
them, captures them on the 2 µm grid and sums them to 8 µm bins, keeping for every bin the cells its counts come from.

PI: Luca Vedovelli. Reviewer: Daniele Sabbatini.

## Status (2026-11-12)

| part | state |
|---|---|
| **Step 1, cell layer** (`R/04b*.R`) | complete and consolidated; six tissue presets calibrated on public 10x Visium HD datasets with H&E; merged into `dev` |
| **Step 2, Visium HD grid** (`R/05b*.R`) | complete and consolidated: counts on 2 µm bins and 8 µm bins with their ground truth, compared with the real reference tissues; whole capture areas simulated in pieces; merged into `dev` |
| next | first analysis methods on simulated and real data (plan in `docs/ROADMAP.md`, section 7) |
| 2025 expression pipeline (`R/01*.R`-`R/03*.R`, `R/04_clustering.R`, `R/05_grid_sampling.R`, `R/06*.R`, `R/test functions/`) | legacy code from before Step 1; not validated against the cell layer and superseded by Step 2. Its README is kept in `docs/legacy/README_2025.md` |

Every biological value used by the presets has a row in `docs/BIO_REFERENCES.md` with a DOI or with the script and the dataset that
estimated it. Nothing is set "from memory".

## The cell layer

`place_cells()` (`R/04b_cell_placement.R`) runs five stages:

1. `extract_regions()` (`R/04b1_extract_regions.R`): pixel clusters become tissue-region polygons in µm.
2. `seed_centroids()` (`R/04b2_seed_centroids.R`): cell centres are drawn in each region at the target density, with a minimum distance
   between centres (uniform seeding by default; optional seeding by cell-type compartments, with voids and density heterogeneity).
3. `tessellate_voronoi()` (`R/04b3_tessellate_voronoi.R`): Voronoi territories of the centres, clipped to the regions.
4. `derive_cell_geometry()` (`R/04b4_derive_cell_geometry.R`): one nucleus per territory, drawn from a size model fitted on real
   nuclei and clipped inside the territory while keeping the target area, up to a cap on the nucleus-to-territory ratio of each preset.
5. `validate_cell_layer()` (`R/04b5_cell_layer_validation.R`): internal checks (geometry validity, area conservation, containment,
   densities and sizes against the preset); with `strict = TRUE` a failed check stops the run.

The six presets live in `R/04b0_tissue_presets.R`. The scale comes from the user, as a number or read from Space Ranger's
`scalefactors_json` (see below). `estimate_pixel_size()` (`R/04b6_estimate_pixel_size.R`) estimates it from the image, but it rejected the image or gave a wrong
scale in five tissues out of six, so `place_cells()` accepts the estimate only with `allow_estimate = TRUE`. The design is `docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md` (version 1.9).

## Tissue presets

One preset per tissue archetype of the test panel. Each preset is calibrated on one public 10x Genomics Visium HD dataset with H&E
(licence CC BY 4.0; map and checksums in `data/real/README.md`). The density is the manual count of nuclei in annotated windows of the
real H&E (row of `docs/BIO_REFERENCES.md` in brackets).

| preset | tissue | reference dataset | target density (nuclei/mm²) |
|---|---|---|---|
| `epithelial` | mouse small intestine | `Visium_HD_Mouse_Small_Intestine` | 12 419 (B-043) |
| `tumor_microenv` | human colorectal cancer | `Visium_HD_Human_Colon_Cancer` | 8 924 (B-044) |
| `stromal_rich` | stroma of the colorectal cancer | `Visium_HD_Human_Colon_Cancer` (stromal regions) | 3 100 (B-045) |
| `lymphoid` | human lymph node | `Visium_HD_Human_Lymph_Node_FFPE` | 28 096 (B-046) |
| `brain_cortex` | mouse brain cortex | `Visium_HD_6p5mm_Mouse_Brain` | 1 185 (B-047) |
| `columnar_muscle` | human heart | `Visium_HD_6p5mm_Human_Heart` | 988 (B-048) |

By default each preset has several cell types (for example epithelial, immune and stromal cells in the intestine), seeded without looking at
the type; the cell-type labels are in Italian (`epiteliale`, `immunitario`, ...). The heart is a planned negative control: isotropic
Voronoi territories cannot reproduce elongated cardiomyocytes.

## The Visium HD grid (Step 2)

`simulate_visium_hd_pieces()` (`R/05b6_simulate_visium_hd_pieces.R`) turns a cell layer into Visium HD data. For each cell it draws a
subtype, a number of transcripts and their genes from expression profiles of the same tissue; it places the transcripts in the nucleus
or in the cytoplasm and moves them (a narrow blur, wider near the tissue edge, an exponential tail of tens of µm for a few per cent of
transcripts, and a second, longer tail that reaches the glass); it
captures them on the 2 µm grid with row and column stripes and a background, and sums them to 8 µm bins. The result holds the counts
(`counts_2um`, `counts_8um`), the bins, and the truth: which cell each captured transcript comes from (`truth$origin`), how much of each
bin each territory and nucleus covers (`truth$overlap`), the background, the transcripts generated for each cell and the technical layers
used.

The area is simulated in square pieces, each cell once, with random draws tied to fixed 56 µm tiles: the same seed gives the same counts
with any piece size and any number of processes, from a small square to a whole capture area. The usage guide,
[`vignettes/visium_hd_simulation.qmd`](vignettes/visium_hd_simulation.qmd) (rendered:
[`vignettes/visium_hd_simulation.html`](vignettes/visium_hd_simulation.html)), covers piece size, time and memory, technical layers,
cropped areas, resuming long runs, the ground truth and the declared limits. The design is
`docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md` (version 0.6).

The expression profiles of the six tissues are read from the lab server's data volume
(`/mnt/micron/geo_spatialtrans/S2.2/profiles`): Step 2 runs on the lab server only, like the test battery.

## Installation

The project is developed and tested on Ubuntu 24.04 (noble) with R 4.6 (4.6.2 on the lab server since 2026-11-05, 4.6.1 before) and the system libraries GDAL 3.8.4, GEOS 3.12.1 and PROJ 9.4.0.
The code is loaded with `source()`: the package is not installed (there is no `NAMESPACE` yet).

```bash
git clone https://github.com/lucavd/2025.geo_spatialtrans.git
cd 2025.geo_spatialtrans
git checkout dev
Rscript --vanilla tools/restore_library.R
```

`tools/restore_library.R` rebuilds the project library `renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu` from `renv.lock`
(199 packages) using Posit Package Manager binaries for Ubuntu noble, then checks every version against the lockfile. Run scripts with
`Rscript --vanilla` and set `.libPaths()` to that folder, as in the quick start below. Other operating systems have not been tried.

## Quick start

A 300 µm × 300 µm square of tissue with the intestine preset. From the repository root, after the installation:

<!-- quickstart:start -->
```r
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)

print(tissue_presets())                       # the six tissue presets

# A 300 x 300 um square of tissue with 1 um pixels, all in one cluster.
clust <- expand.grid(x = 1:300, y = 1:300)
clust$intensity_cluster <- 1L

cl <- place_cells(clust, tissue_preset = "epithelial", pixel_size_um = 1, random_seed = 42)

nrow(cl$cell_df)                              # number of simulated cells
head(cl$cell_df)
cl$region_check                               # area of the region and of its territories
p <- plot_cell_layer(cl, mode = "voronoi")
ggplot2::ggsave("cell_layer_epithelial.png", p, width = 6, height = 6, dpi = 150)
```
<!-- quickstart:end -->

`place_cells()` prints one line such as `[place_cells] epithelial: 1 regioni, 1118 cellule, scala 1 um/px (user), controlli 0 ERROR / 0 WARN, 1.8 s`.
The same seed gives the same cells, bit for bit.

### Quick start, Step 2

The same kind of square, simulated as Visium HD data. The square is cut out of a tissue that continues (`frame_open = TRUE`); on an area
this small a single piece (`piece_um = Inf`) is the fastest choice. It takes about half a minute on the lab server.

<!-- quickstart2:start -->
```r
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(Matrix) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)

clust <- expand.grid(x = 1:150, y = 1:150)    # 300 x 300 um of tissue, 2 um pixels
clust$intensity_cluster <- 1L
cl <- place_cells(clust, tissue_preset = "epithelial", pixel_size_um = 2, random_seed = 42)

hd <- simulate_visium_hd_pieces(cl, seed = 1, piece_um = Inf, frame_open = TRUE)
dim(hd$counts_2um)                            # genes x 2 um bins (tissue and 280 um of glass around it)
dim(hd$counts_8um)                            # genes x 8 um bins
summary(Matrix::colSums(hd$counts_8um)[hd$bins_8um$in_tissue])   # UMI per 8 um bin in the tissue
head(hd$truth$origin)                         # captured transcripts per (2 um bin, cell)
all(unlist(hd$metadata$checks))               # internal checks of the run
```
<!-- quickstart2:end -->

### Input from an image

`clust` is a data frame of tissue pixels: integer `x` and `y` and a cluster column (`intensity_cluster` by default; set `cluster_col`).
`pixel_size_um` is required: a number in µm per pixel, or `list(scalefactors_json = <file>)` to read it from Space Ranger. Each cluster
becomes one or more tissue regions. The 2025 functions `generate_synthetic_tissue()` (`R/03b_generate_synthetic_tissue.R`) and
`cluster_image()` (`R/04_clustering.R`) produce this format and built the synthetic inputs of the first cell-layer tests
(`tools/S1.1_inputs.R`).

### Main options of `place_cells()`

| argument | default | meaning |
|---|---|---|
| `tissue_preset` | `"epithelial"` | one of `tissue_presets()` |
| `random_seed` | `42` | seed of the three stages (seeding, nucleus size, pale nuclei; `stage_seeds()`) |
| `types` | `"multi"` | several cell types per preset; `"single"` = one type, as in the first calibration |
| `spatial` | preset default (`"uniform"`) | `"compartments"` seeds cell types in spatial compartments, with voids |
| `cell_types`, `region_composition` | `NULL` | a user catalogue of cell types and the composition of each cluster, instead of a preset |
| `validate`, `strict` | `TRUE`, `TRUE` | run the internal checks; stop on a failed check |

The full list, with the options kept to reproduce earlier steps, is the signature of `place_cells()` in `R/04b_cell_placement.R`; its roxygen
header documents most arguments but not yet all of them (for example `random_seed`, `validate` and `strict`).

### Output

A `cell_layer` list (design §4.1): `cell_df` (one row per cell: position, type, territory area, nucleus geometry and area),
`cell_territories` and `cell_nuclei` (`sf` polygons), `region_df` and `region_polygons`, `territory_df`, `region_check`,
`region_composition`, `cell_types_resolved` and `metadata` (seeds, scale, check results, run time). `plot_cell_layer()` draws territories,
types, density or regions.

## Known limits of the simulated cells

Declared in the design (§7.4) after comparison with the real tissues; keep them in mind when testing methods on the simulated data:

- **too many nuclei touching each other**: shown in the stroma of the colorectal cancer and for pairs of nuclei, undecided in the
  colorectal cancer, no longer shown in the lymph node (`docs/decisions/DEC-063.md`);
- **territories too uniform** with the default uniform seeding (territory-area variability and number of neighbours below the real ones in
  the intestine, the colorectal cancer and the lymph node); seeding by compartments with voids corrects part of it but is not the default;
- **heart**: isotropic territories, while real cardiomyocytes are elongated (planned negative control);
- **preset densities** refer to tissue fully covered by cells, as in the manual count, not to tissue with its voids.

The limits of the simulated counts (variability outside the nucleus, far glass of the colorectal cancer, empty bins of the dense tissues,
fibre orientation of the heart, time of the largest sections) are listed in the usage guide, section «Known limits», and in section 17
of the Step 2 design.

## Tests

The test battery runs on the lab server only. Besides the code it reads the real datasets on the server's data volume, derived files of
the reference steps that are too large for git (tissue masks and nuclei of the real regions, annotation windows) and the project's Python
environment. In a new clone on the server, after `tools/restore_library.R`, link these inputs from the server's working copy:

```bash
bash tools/link_server_inputs.sh      # symbolic links to every git-ignored input of the tests (libraries, logs and caches excluded)
bash R/testing/run_S1.6.sh test       # full regression battery of Step 1 (PASS/FAIL per assertion, about 35 minutes)
bash R/testing/run_S1.6.sh mutants    # altered copies of the code that the tests must catch
bash R/testing/run_S1.6.sh env        # DESCRIPTION vs code, restore from the lockfile into an empty library
for s in 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8; do Rscript --vanilla R/testing/test_S$s.R; done   # Step 2 batteries (PASS/FAIL per assertion)
```

The Step 1 battery prints `GATE G-1 PASS` when every assertion passes except four known failures recorded as pre-registration errors
(`docs/decisions/DEC-011.md`, `docs/decisions/DEC-023.md`). The Step 2 batteries print `PASS` or `FAIL` for each assertion; they
take from a few minutes to about 25 minutes each and can run in parallel. Count the `FAIL` lines: `R/testing/test_S2.3.R` exits with
status 0 even when an assertion fails. The Step 1 quick start needs neither the data nor the links; the
Step 2 quick start and the usage guide need the expression profiles on the data volume.

## Repository map

| path | content |
|---|---|
| `R/04b*.R` | the cell layer (Step 1) |
| `R/05b*.R` | the Visium HD grid (Step 2): grid and geometric truth (`R/05b1_build_hd_grid.R`), profiles (`R/05b2_simulate_cell_profiles.R`), position and movement of transcripts (`R/05b3_place_and_move_transcripts.R`), capture and binning (`R/05b4_capture_and_bin.R`), the reference engine (`R/05b5_simulate_visium_hd.R`) and the entry point (`R/05b6_simulate_visium_hd_pieces.R`) |
| `vignettes/` | usage guide of Step 2 (Quarto source and rendered HTML) |
| `R/testing/` | tests (`R/testing/test_<step>.R`), gate scripts (`R/testing/run_<step>.sh`), mutants (`R/testing/mutants_<step>.py`) |
| `tools/` | analysis, calibration and benchmark scripts of each step; `tools/okf.py` checks the documentation bundle |
| `results/<step>/` | versioned results tables of each step |
| `reports/<step>.qmd` | step reports in Italian (Quarto source; the rendered `reports/<step>.html` sits next to it) |
| `docs/superpowers/specs/` | design documents of Step 1 and Step 2 |
| `docs/BIO_REFERENCES.md` | biological references with sources |
| `docs/ROADMAP.md`, `docs/BACKLOG.md`, `CHANGELOG.md` | plan, open items, history |
| `docs/steps/`, `docs/decisions/`, `docs/index.md` | one file per step and per decision (Open Knowledge Format), index |
| `data/real/README.md` | the real reference datasets, their checksums and scale |

Project documents are in Italian; code, tables and figures are in English.

## Licence and authors

MIT licence (`LICENSE`, copyright Luca Vedovelli). Authors and roles are in `DESCRIPTION`.
