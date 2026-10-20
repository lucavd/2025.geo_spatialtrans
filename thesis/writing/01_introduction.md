# Introduction

<!--
Skeleton of bullets to expand later. Target: 8–10 pages (PROTOCOLLO.md §8); first 2–3 page draft due week 1.
Order follows §8: spatial transcriptomics and deconvolution → simulation in benchmarks → Crowell 2023 → generative models → the question.
[REF] = citation still to add (no claims or numbers without a source). Limits go in the Discussion, not here.
-->

## 1. Spatial transcriptomics and the spot

- Spatial transcriptomics measures gene expression while keeping track of *where* in the tissue it was measured. [REF]
- In spot-based technologies (e.g. Visium), each measurement comes from a **spot**: a small capture area that collects the RNA of several cells at once. [REF]
- So a spot is not a cell. It is a **mixture**: its counts are the sum of the transcripts of the cells that fall in it.
- Two dimensions of "spatial" are worth separating:
  - **geometric:** positions, neighbouring spots, tissue structure;
  - **spot (mixture):** each measurement is a blend of cell types in unknown proportions.
- This thesis works on the **spot dimension**: how well can we recover what is inside each spot?

## 2. Deconvolution: unmixing the spot

- **Deconvolution** estimates the proportion of each cell type in each spot. [REF]
- Most methods need a **reference**: the typical expression profile of each cell type, usually taken from single-cell RNA-seq.
- Methods used here: RCTD, CARD, cell2location, NNLS as a simple baseline (SPOTlight optional). [REF for each]
  - one line each: statistical model, what reference it needs.
- The quality of the result depends on the method *and* on how well the reference matches the cells in the spots.

## 3. Why benchmarks need simulation

- To score a deconvolution method you need the **true composition** of each spot, which real tissue does not give you.
- Simulation solves this: build spots with a known recipe, then check how close each method gets. [REF]
- Common designs: pseudo-spots made by summing single cells; fully simulated counts. [REF]
- Two families of simulators: **de novo** (counts generated from parameters the user chooses: flexible, but possibly unrealistic) and **reference-based** (parameters estimated from a real dataset: more realistic, but limited to what the reference contains). Many are hybrids: a baseline fitted to one homogeneous group, with differences between groups added by hand. [Crowell 2023]
- Simulation studies should be planned like experiments: aims, data-generating mechanism, estimands, methods, performance measures (ADEMP, Morris 2019). [REF]
- **The catch:** a benchmark's ranking can depend on choices in the simulator, not only on the methods.

## 4. Crowell 2023: simulators are not neutral

- Crowell et al. compared 16 scRNA-seq simulators that learn from real data (reference-based) on 12 real datasets, then re-ran two benchmarks (batch correction, clustering) on real and simulated versions of the same data. [Crowell 2023; Author Correction 2024 fixes only a table legend]
- **Fidelity is limited.** Most simulators can generate only one group of cells unless the user adds artificial differences; the structure *between* groups (batches, clusters) is the worst-reproduced feature. [Crowell 2023]
- **Benchmarks become over-optimistic.** Clustering methods scored higher on simulated data than on the real data the simulation was fitted to, for every simulator. [Crowell 2023]
- **Rankings move.** Clustering rankings on simulated data agreed poorly with those on real data; among simulators of similar quality, they depended more on which reference dataset was used than on which simulator. For batch correction the results were closer to real data. [Crowell 2023]
- **Realistic is not the same as reliable.** The most faithful simulators did not give the most faithful rankings (e.g. muscat: high fidelity, inconsistent rankings; SPsimSeq: mediocre fidelity, the best integration ranking). Which features of the data matter for a given task is an open question. [Crowell 2023]
- Their study covered single-cell tasks only; spot deconvolution was not tested.
- Link: if the simulator can change conclusions, and overall realism does not tell us which simulator to trust, we need to isolate **one ingredient** and measure its effect on the ranking directly. For deconvolution, the obvious ingredient is the cell-type expression profiles (§5).

## 5. Where cell-type expression comes from

- The one ingredient we vary: the **cell-type expression profiles** that go into the spots.
- Possible sources:
  - **empirical:** averages from real single-cell data;
  - **parametric:** profiles drawn from a statistical distribution, with controlled marker genes;
  - **simulator-based:** e.g. scDesign3 [REF];
  - **generative AI models:** e.g. GEM-1, which produces expression from metadata alone. [REF: GEM-1 preprint]
- Generative models are new, and their training data and evaluation have limits (see `letture.md`).

## 6. The question

- **Main question:** does the source of cell-type expression change the rankings of a deconvolution benchmark?
- Objective 1 (fidelity): how close is each source to real held-out cells?
- Objective 2 (benchmark sensitivity): do method rankings and errors change across sources, and which source gives the ranking closest to pseudo-spots built from real held-out cells?
- The two objectives are measured separately on purpose: Crowell 2023 showed that a more realistic simulator does not necessarily give a more faithful ranking, so Objective 1 cannot stand in for Objective 2.
- Controls: a negative control that must fail (permuted genes) and a positive control as the ceiling (a second real donor).
- Three tissues: colon/colorectal cancer, lymphoid tissue, brain. More than one, because benchmark rankings can depend on the reference dataset as much as on the simulator. [Crowell 2023]
- One sentence on why it matters: guidance for anyone who builds a deconvolution benchmark.
