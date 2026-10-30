# geo.spatialtrans: a Visium HD simulator with a cell layer

`geo.spatialtrans` simulates Visium HD spatial transcriptomics data in order to test analysis methods against a known truth.
Its distinctive feature is an explicit **cell layer** between the tissue image and the 2 µm capture grid: every simulated cell has a
position, a territory (a Voronoi tile clipped to its tissue region) and a nucleus, so that later steps can map transcripts from cells to
2 µm spots and 8 µm bins with a known ground truth.

PI: Luca Vedovelli. Reviewer: Daniele Sabbatini.

## Status (2026-10-30)

| part | state |
|---|---|
| **Step 1, cell layer** (`R/04b*.R`) | complete and consolidated; six tissue presets calibrated on public 10x Visium HD datasets with H&E; merged into `dev` |
| **Step 2, Visium HD grid** (expression from cells to 2 µm spots and 8 µm bins) | not started; design in preparation |
| 2025 expression pipeline (`R/01`-`R/03`, `R/04_clustering.R`, `R/05`-`R/06k`, `R/test functions/`) | legacy code from before Step 1; not validated against the cell layer; it will be replaced by Step 2. Its README is kept in `docs/legacy/README_2025.md` |

Every biological value used by the presets has a row in `docs/BIO_REFERENCES.md` with a DOI or with the script and the dataset that
estimated it. Nothing is set "from memory".

## The cell layer

`place_cells()` (`R/04b_cell_placement.R`) runs five stages:

1. `extract_regions()` (`R/04b1_extract_regions.R`): pixel clusters become tissue-region polygons in µm.
2. `seed_centroids()` (`R/04b2_seed_centroids.R`): cell centres are drawn in each region at the target density, with a minimum distance
   between centres (uniform seeding by default; optional seeding by cell-type compartments, with voids and density heterogeneity).
3. `tessellate_voronoi()` (`R/04b3_tessellate_voronoi.R`): Voronoi territories of the centres, clipped to the regions.
4. `derive_cell_geometry()` (`R/04b4_derive_cell_geometry.R`): one nucleus per territory, drawn from a size model fitted on real
   nuclei and clipped inside the territory while keeping the target area.
5. `validate_cell_layer()` (`R/04b5_cell_layer_validation.R`): internal checks (geometry validity, area conservation, containment,
   densities and sizes against the preset); with `strict = TRUE` a failed check stops the run.

The six presets live in `R/04b0_tissue_presets.R`; `estimate_pixel_size()` (`R/04b6_estimate_pixel_size.R`) reads the image scale when the
user does not give it. The design is `docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md` (version 1.8).

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

## Installation

The project is developed and tested on Ubuntu 24.04 (noble) with R 4.6.1 and the system libraries GDAL 3.8.4, GEOS 3.12.1 and PROJ 9.4.0.
The code is loaded with `source()`: the package is not installed (there is no `NAMESPACE` yet).

```bash
git clone https://github.com/lucavd/2025.geo_spatialtrans.git
cd 2025.geo_spatialtrans
git checkout dev
Rscript --vanilla tools/restore_library.R
```

`tools/restore_library.R` rebuilds the project library `renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu` from `renv.lock`
(199 packages) using Posit Package Manager binaries for Ubuntu noble, then checks every version against the lockfile. Every script of the
project calls `.libPaths()` on that folder and runs with `Rscript --vanilla`. On other systems the restore falls back to source packages.

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

`place_cells()` prints one line such as `[place_cells] epithelial: 1 regioni, 1118 cellule, scala 1 um/px (user), controlli 0 ERROR / 0 WARN`.
The same seed gives the same cells, bit for bit.

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

The full list, with the options kept to reproduce earlier steps, is in the roxygen header of `place_cells()`.

### Output

A `cell_layer` list (design §4.1): `cell_df` (one row per cell: position, type, territory area, nucleus geometry and area),
`cell_territories` and `cell_nuclei` (`sf` polygons), `region_df` and `region_polygons`, `territory_df`, `region_check`,
`region_composition`, `cell_types_resolved` and `metadata` (seeds, scale, check results, run time). `plot_cell_layer()` draws territories,
types, density or regions.

## Known limits of the simulated cells

Declared in the design (§7.4) after comparison with the real tissues; keep them in mind when testing methods on the simulated data:

- **too many nuclei touching each other** in dense tissues: colorectal cancer and lymph node, borderline in the stroma of the cancer;
  none of the nucleus models tried brings them to the real values;
- **territories too uniform** with the default uniform seeding (territory-area variability and number of neighbours below the real ones in
  the intestine, the colorectal cancer and the lymph node); seeding by compartments with voids corrects part of it but is not the default;
- **heart**: isotropic territories, while real cardiomyocytes are elongated (planned negative control);
- **preset densities** refer to tissue fully covered by cells, as in the manual count, not to tissue with its voids.

## Tests

The test battery runs on the lab server only. Besides the code it reads the real datasets on the server's data volume, derived files of
the reference steps that are too large for git (tissue masks and nuclei of the real regions, annotation windows) and the project's Python
environment. In a new clone on the server, after `tools/restore_library.R`, link these inputs from the server's working copy:

```bash
bash tools/link_server_inputs.sh      # symbolic links to every git-ignored input of the tests (libraries, logs and caches excluded)
bash R/testing/run_S1.6.sh test       # full regression battery of Step 1 (PASS/FAIL per assertion, about 20 minutes)
bash R/testing/run_S1.6.sh mutants    # altered copies of the code that the tests must catch
bash R/testing/run_S1.6.sh env        # DESCRIPTION vs code, restore from the lockfile into an empty library
```

The battery prints `GATE G-1 PASS` when every assertion passes except four known failures recorded as pre-registration errors
(`docs/decisions/DEC-011.md`, `docs/decisions/DEC-023.md`). The quick start above needs neither the data nor the links.

## Repository map

| path | content |
|---|---|
| `R/04b*.R` | the cell layer (Step 1) |
| `R/testing/` | tests (`test_<step>.R`), gate scripts (`run_<step>.sh`), mutants (`mutants_<step>.py`) |
| `tools/` | analysis, calibration and benchmark scripts of each step; `okf.py` checks the documentation bundle |
| `results/<step>/` | versioned results tables of each step |
| `reports/<step>.qmd`, `.html` | step reports (Italian) |
| `docs/superpowers/specs/` | design of Step 1 |
| `docs/BIO_REFERENCES.md` | biological references with sources |
| `docs/ROADMAP.md`, `docs/BACKLOG.md`, `CHANGELOG.md` | plan, open items, history |
| `docs/steps/`, `docs/decisions/`, `docs/index.md` | one file per step and per decision (Open Knowledge Format), index |
| `data/real/README.md` | the real reference datasets, their checksums and scale |

Project documents are in Italian; code, tables and figures are in English.

## Licence and authors

Authors and roles are in `DESCRIPTION`. The licence declared there is MIT.
