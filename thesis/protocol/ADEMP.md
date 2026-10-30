# ADEMP protocol — expression sources and deconvolution benchmarks

**Status:** DRAFT v0.4, 2026-10-30 (v0.1 of 2026-10-26 with Luca's answers; v0.2 of 2026-10-27 donor counts and noise ceiling; v0.3 of 2026-10-29 split, cell types, E′, noise model, pilot scope, packages; v0.4 Luca's answers of 2026-10-30: cutoff criterion, entry rule for G, held-out group of 4 donors, donor covariate in φ, the `[LUCA]` items; factor levels). Not approved, not pre-registered. Becomes valid only when Luca approves it with a date (Week 2 checkpoint).
**Sources:** `PROTOCOLLO.md` §1, §2, §5; `DECISIONS.md`; Morris, White & Crowther 2019 (doi:10.1002/sim.8086), notes in `writing/letture.md`.
**Placeholders:** `[PILOT]` = fixed from the Week 2 pilot · `[REAL DATA]` = fixed from the training split of the real data · `[LUCA]` = decision for Luca.
Any change after approval is a dated amendment (section 12), written before the results it could affect.

---

## 1. Aims

- **Aim 1 — Fidelity.** How close are the cell-type mean profiles μ (genes × cell types) of each expression source to real **held-out** profiles of the same tissue and cell type?
- **Aim 2 — Benchmark sensitivity.** Do the error and the ranking of deconvolution methods change with the expression source? Which source gives the ranking closest to the one obtained on pseudo-spots built from **real held-out cells** (the real anchor)?

### Hypotheses
- **H1 (negative control).** Arm C− has the lowest fidelity in every tissue, and its method ranking differs from the anchor.
- **H2 (donor-to-donor reference).** Arm E′ shows how close a second real donor gets to the held-out donor, in fidelity and in ranking agreement with the anchor. It gives the scale on which the other arms are read: a source close to E′ is "as good as another real donor". It is not a ceiling. A synthetic arm may reach or exceed E′ without anything being wrong: M (fitted on the training donors) and G (trained on public archives) are estimated from many donors, so they give a population-average profile, which can be closer to a new donor than one other donor is. Exceeding E′ is therefore not evidence of leakage. The ceiling of Aim 1 is the split-half noise ceiling of the held-out donor (5.1); leakage is checked on the provenance of the data (section 9). *(DECISIONS 2026-10-26: E′ as donor-to-donor reference, Pierangela, approved by Luca; noise ceiling and provenance-based leakage check, Luca.)*
- **H3 (main question).** At least one synthetic arm (P, M, G) gives a method ranking that differs from the anchor ranking.
- "Differs" and "closest" are turned into decisions in section 5.4.

## 2. Data-generating mechanisms

### 2.1 Main factor: source of μ
| Arm | Source of μ | Role |
|---|---|---|
| E | Empirical means from the real scRNA training split | "Honest" reference |
| P | Parametric de novo: means from a gamma distribution, markers with controlled fold-changes | Maximum control, minimum realism |
| M | scDesign3 trained on the training split | Model-based middle ground |
| G | GEM-1-sc, "mean estimation" mode, from CL/UBERON metadata | Generative model (optional; enters by the rule of DECISIONS 2026-10-30) |
| C− | E profiles with gene labels permuted within each cell type | Negative control: must fail |
| E′ | Means from a second, independent real donor ("independent comparator donor"; nothing is tuned on it, so it is not a validation set) | Reference for donor-to-donor variability (not a ceiling). In Aim 2 its means generate the spots like every arm; the method reference is always the training scRNA (PROTOCOLLO §5.2) (DECISIONS 2026-10-29, approved by Luca 2026-10-30) |

- P parameters (gamma shape/rate, number of markers, fold-change): `[REAL DATA]` (matched to the E arm's mean–variance range).
- C− permutation: one new permutation per replicate, so the result does not depend on one particular shuffle (DECISIONS 2026-10-30, Luca).
- G is evaluated in all three tissues. Criterion (DECISIONS 2026-10-30, Luca, replacing the unapproved 2026-10-27 row): GEM-1-sc was trained on the whole CELLxGENE of July 2024 (preprint, doi:10.1101/2025.09.08.674753, Methods), so a dataset is usable for G if it is **not in the CELLxGENE Census of 2024-07-01**. Luca checked that Census version: PsychAD, Pelka 2021 and Wells 2025 are all absent. Potential limit (§11): the CELLxGENE version behind the gem-1-sc served today by the API is not known.
- **Entry of G** (DECISIONS 2026-10-30, Luca): no deadline at 30/10 any more. G enters when written permission from Synthesize Bio, API access and vocabulary coverage are all in place; the decision is written in DECISIONS, with the date, before any G result is looked at. If G cannot arrive in time for the thesis, it is removed and declared as a limit. Cut rule if late: arm G, then the brain tissue (DECISIONS 2026-10-26).

### 2.2 Secondary factors (fully crossed, fixed grid)
| Factor | Levels |
|---|---|
| Tissue | colon/CRC, lymphoid, brain |
| Cells per spot | 5 and 15 (Pierangela, 2026-10-30): 5 = middle of the 2–8-cell synthetic mixtures of SPOTlight (Elosua-Bayes et al. 2021, doi:10.1093/nar/gkab043); 15 = middle of the 10–20 cells of a 50–100 µm spot (same source). A fixed number per level, so the methods that accept it get the true value (oracle information, DECISIONS 2026-10-30) |
| Sequencing depth | 10,000 and 20,000 UMI per spot (Pierangela, 2026-10-30): 10,000 = RCTD's simulated pixels "to model the conditions of Visium" (Cable et al. 2022, doi:10.1038/s41587-021-00830-w); 20,000 = SPOTlight and synthspot (Spotless), which downsample synthetic spots around 20,000 UMI. Every arm and the anchor are thinned at random to the target; a spot below the target is kept as it is (synthspot rule). `03_reference.R` reports how often 5 or 15 cells fall below each target |

Full design: 5 arms (E, P, M, C−, E′) × 3 tissues × 2 × 2 = **60 DGMs**; **72** if arm G enters (6 arms, all three tissues; §2.1). Plus the real anchor for each tissue × cells-per-spot × depth (12 conditions).
- Number of cell types per tissue and their labels: §2.6 (from `02_split.R`).
- Number of spots per pseudo-dataset: `[PILOT]`. The held-out donors are typical, hence small (e.g. lymphoid 583B, 2,433 cells; colon C104, 2,886), so pseudo-spots and the anchor reuse the same held-out cells many times across spots and replicates, as in other pseudo-spot simulators (e.g. synthspot in Spotless). The pilot reports the expected reuse per cell for the chosen number of spots, and the choice is written here.
- Composition (Pierangela, 2026-10-30): π_s ~ symmetric Dirichlet(α = 1) over the spot types, then the cells of the spot drawn multinomially; the truth is the cell fraction. α = 1 is the default of the Open Problems spatial-decomposition simulator (GitHub `openproblems-bio/task_spatial_decomposition`, `dataset_simulator`, which also offers 0.5 and 5). synthspot (Spotless) uses tissue-region priors instead, which pseudo-spots without regions do not have.
- Spots per pseudo-dataset in the pilot: 500 (synthspot `generate_synthetic_visium_lite` default); final number `[PILOT]`.

### 2.3 Generating one replicate
1. For each spot, draw a composition from a Dirichlet(α) and sample the cells.
2. For each cell of type k, the count of gene g is negative binomial with mean μ_gk (scaled by library size) and dispersion φ_g estimated from the real training data. **The noise model is the same for all arms**; only μ changes.
   Estimation (DECISIONS 2026-10-29; donor covariate in φ by Luca, 2026-10-30; sources in `notes/noise_model_digest.md`):
   - μ_gk = relative mean profile of type k on the training split, pooled: sum of the counts of gene g over the cells of type k / sum of all their counts (each column sums to 1). This is the maximum-likelihood rate of a Poisson model with log library size as offset. Same "profile × library size" structure as RCTD (Cable et al. 2022, doi:10.1038/s41587-021-00830-w: "average normalized expression of gene j within all cells of cell type k", counts ~ Poisson(N_i λ)) and muscat (Crowell et al. 2020, doi:10.1038/s41467-020-19894-4: mean = relative abundance × library size). RCTD averages per-cell ratios instead; the pooled version weights cells by library size and matches our own noise model.
   - Library size L of a simulated cell: drawn at random from the real library sizes of training cells of the same type; mean count = μ_gk × L. Empirical, as in muscat (simulated cells take the library size of reference cells of their group), rather than parametric (log-normal in Splatter, Zappia et al. 2017, doi:10.1186/s13059-017-1305-0).
   - φ_g: one per gene, fitted on the training split with `glmGamPoi` (design `~ cell_type + donor`, L as offset; the donor covariate keeps donor-to-donor variability out of φ, otherwise the simulated arms would be noisier than the anchor, which comes from one donor: DECISIONS 2026-10-30, Luca; Ahlmann-Eltze & Huber 2020, doi:10.1093/bioinformatics/btaa1009: better overdispersion estimates "on datasets with many small counts", no dataset size limit), on a sample of 20,000 training cells in every tissue, stratified by cell type (proportional to each type's share) with a fixed seed: one rule for all tissues, because all brain cells (~900,000) do not fit in memory; `glm_gp()` itself offers subsampling (`subsample`). Pre-registered check in colon, where all training cells fit: φ from the sample vs φ from all cells; the sample is adequate if the correlation of log φ across genes is ≥ 0.99 and the median ratio φ_sample/φ_all is between 0.95 and 1.05, otherwise the sample rises to 50,000 cells in every tissue (DECISIONS 2026-10-29, refined 2026-10-30; `tests/test_phi_sample.R`). Per-gene NB dispersion as in muscat (edgeR `estimateDisp`) and in cell2location's model (Kleshchevnikov et al. 2022, doi:10.1038/s41587-021-01139-4: "gene- and batch-specific over-dispersion").
3. Sum the cells of the spot.

**Real anchor:** pseudo-spots built by summing **real cells** of the held-out donor, with the same Dirichlet compositions. The truth is known because each cell's type is known.
- Every anchor pseudo-spot is built inside a single sample, because a real spot never mixes cells from two samples (DECISIONS 2026-10-26). Colon uses MMRp tumour samples only and lymphoid the thoracic lymph node only (DECISIONS 2026-10-26 "Dataset", refined 2026-10-27, approved by Luca 2026-10-30), so there is one μ per tissue.
- Depth of the anchor: its spots are thinned at random to the same target depth as the simulated arms (DECISIONS 2026-10-30, Luca).

### 2.4 Common random numbers
This follows from PROTOCOLLO §5.2 (only the source of μ changes; the anchor uses "la stessa composizione Dirichlet"). Within replicate i, all arms and the anchor use the **same** compositions and the same cell-type draws per spot; only μ (and therefore the counts) differs. Differences between arms are then not inflated by different compositions, which lowers the Monte Carlo SE of between-arm comparisons (Morris §5.3). All methods are run on the same datasets (Morris recommendation).

### 2.5 Data split
Split by donor, never by cell (DECISIONS 2026-10-14 and 2026-10-26). One dataset per tissue; training, held-out and E′ are different donors of the same dataset (DECISIONS 2026-10-26). Training split → E, M, φ_g, P calibration and the reference given to every method. Held-out donor → Aim 1 target, noise ceiling and the anchor. Second independent donor → E′.
- Colon: MMRp tumour samples only (DECISIONS 2026-10-26 "Dataset", refined 2026-10-27, approved by Luca 2026-10-30). All samples of a donor go to the same split (e.g. C130 has two tumour lesions, C130_TA and C130_TB): no pseudoreplication, no leakage.
- Lymphoid: thoracic lymph node only (DECISIONS 2026-10-26 "Dataset", refined 2026-10-27, approved by Luca 2026-10-30).

| Tissue | Dataset | Chemistry | Donors: training | held-out | E′ |
|---|---|---|---|---|---|
| Brain | PsychAD (DLPFC), healthy adults: disease "normal", age ≥18 years (208 donors; 201 after donor QC, `02_split.R`) | 10x 3′ v3, single-nucleus | 197 (all others) | 1 (+2 in the Aim 1 group) | 1 |
| Colon/CRC | Pelka 2021, 3′ v2 subset, MMRp tumours (20 donors) | 10x 3′ v2 | 16 | 1 (+2) | 1 |
| Lymphoid | Wells 2025, thoracic lymph node (21 donors) | 10x 5′ v2 | 17 | 1 (+2) | 1 |

- One held-out donor and one E′ donor per tissue, as in PROTOCOLLO §5.2 ("il donatore held-out", "un secondo donatore reale indipendente"). Same rule in every tissue: after the biological selection (healthy adults; MMRp tumours; thoracic lymph node), all other eligible donors go to training; no donor is dropped at random. Training sizes therefore differ between tissues, which does not affect the analyses, made separately by tissue (declared limit, §11). Counts from the CELLxGENE metadata (obs) of the three datasets, 2026-10-27; `02_split.R` recomputes them after the §2.6 filters.
- **Donor QC** (DECISIONS 2026-10-29, approved by Luca 2026-10-30): a donor is excluded from the whole dataset when its number of cells is a low outlier, below median − 3 MAD on the log scale (the `scater::isOutlier` convention for cells, OSCA, Amezquita et al. 2020, doi:10.1038/s41592-019-0654-x, applied to donors). Excluded: brain 7 donors (61–770 nuclei, 3,335 in all; lower bound 876), colon and lymphoid none.
- **Which donors** (rule fixed before any result, DECISIONS 2026-10-29, approved by Luca 2026-10-30): among the donors left after QC, donors are drawn at random with the fixed seed 20261029 (`R/assign_split.R`); the first is held-out, the second E′. A random draw makes the held-out a typical donor rather than the best-sampled one. Result: brain Donor_207 (7,271 nuclei, 20 types with ≥50 cells) / Donor_680 (3,508; 11); colon C104 (2,886; 8) / C105 (3,281; 7); lymphoid 583B (2,433; 9) / 591C (10,775; 13).
- **Held-out group** (DECISIONS 2026-10-30, Luca): 4 donors are drawn with the same seed; the first two are unchanged (checked by `02_split.R`). Aim 1 is computed on each of the 4 donors and E′ on all 6 pairs of the group, giving a distribution of donor-to-donor variability instead of one number. Aim 2 uses the first donor (held-out) with the second as E′; a second held-out donor as sensitivity analysis if the pilot shows there is time. Donors 3–4: brain Donor_307 (5,181; 16) / Donor_957 (3,873; 14); colon C134 (2,563; 7) / C162 (7,450; 13); lymphoid D543 (9,254; 12) / D512 (7,375; 11). None of the three datasets is in CELLxGENE of July 2024 (§2.1), so every donor is usable for G.
- **Compute fallback (brain):** if the pilot shows that a step cannot handle the full training split (~1 million nuclei; e.g. scDesign3 for arm M, the cell2location reference), cells (not donors) are subsampled within the training split with a fixed seed, and this is written as a dated amendment (§12) before any result.

Donor counts per split are filled by `scripts/02_split.R` (Luca's request, 2026-10-26); candidates and rules in `data_inventory.xlsx`, sheet `selected`.

### 2.6 Filtering of the real data (before building μ) — proposal `[LUCA]`
PROTOCOLLO §5.2 does not mention this step; it is proposed to Luca. The same minimal filter is applied to every real dataset before μ, φ_g, the reference for the methods and the anchor are computed. It follows Crowell et al. 2023 (Additional file 1, §6.1.2), who filtered their real references "as insufficient numbers of cells and low-quality observations can interfere with estimation of simulation parameters". Filtered data are saved in `data/derived/`; the downloaded files are never changed (PROTOCOLLO §3.7).
1. **Cell types:** remove a cell type from a split (training, held-out donor, E′ donor) if it has fewer than 50 cells there. In Crowell the groups are batches or clusters; for us they are cell types. The cell types removed are listed in the report.
   Consequences (DECISIONS 2026-10-29, approved by Luca 2026-10-30): sources E, P, M, G, C− and the method reference use every type with ≥50 training cells (brain 27, colon 18, lymphoid 30; `type_kept` in `<tissue>_types.csv`). Aim 1 scores only the types with ≥50 held-out cells (20 / 8 / 9); the others are "not evaluated", not penalised. Pseudo-spots of every arm and the anchor use the types with ≥50 cells in training, held-out and E′ (brain 11, colon 6, lymphoid 7), because E′ has no profile for some held-out types (9 / 2 / 2): every arm, E′ included, then gets the same compositions (§2.4). The method reference keeps all training types, so reference types absent from the spots have true proportion 0 and proportion given to them counts as error (false positive). E′ has 11 / 7 / 13 types; in Aim 1 it is compared with the held-out on the types with ≥50 cells in both (11 / 6 / 7).
2. **Genes:** keep genes with a count greater than 1 in at least 10 cells. Proposal: decide the gene list on the training split only and use the same list for every split and arm, so that all μ matrices have the same genes and the held-out data do not influence the choice.
3. **Cells:** keep cells with at least 100 detected genes (non-zero counts).

## 3. Estimands

- **Primary (Aim 2):** the proportion π_sk of each cell type k in each spot s. It is a vector per spot, not a single scalar θ as in Morris; we score it with RMSE and Jensen–Shannon divergence (section 5).
- **Ranking target (Aim 2):** the ordering of the methods by mean error within each arm × condition, compared with the ordering on the anchor.
- **Secondary (Aim 1):** the cell-type profiles μ and their markers, compared with the held-out real profiles.

## 4. Methods

| Method | Package / language | Notes |
|---|---|---|
| RCTD | spacexr (R; v2.2.1, GitHub dmcable/spacexr commit 9f5dc33, installed 2026-10-29) | mode: **full** (DECISIONS 2026-10-30, Luca). Cable et al. 2022 recommend doublet mode only when pixels rarely hold three or more cell types (Slide-seq), otherwise full mode (e.g. Visium); our pseudo-spots mix several cell types, so `doublet_mode = "full"` (package default is "doublet", so it must be set). No coordinates needed: `SpatialRNA(..., use_fake_coords = TRUE)` (checked in the source, 2026-10-23). See `letture.md` |
| cell2location | Python, GPU on the server | n_sim may be lower; fixed by dated amendment in the first days on the server (DECISIONS 2026-10-14). Needs raw counts. Hyper-parameters: `N_cells_per_location` = the true cells-per-spot level of each condition (§2.2; oracle information, DECISIONS 2026-10-30); `detection_alpha`: proposal (Pierangela, 2026-10-30, to Luca) 200 as primary, 20 as sensitivity. The official tutorial (checked 2026-10-30) does not tie the choice to the technology of reference and spots: 20, the current default, is for data with "large technical variability in RNA detection sensitivity" (typical of human Visium), 200 for datasets with "low technical effects", and it recommends testing both; pseudo-spots have no spatial technical effect |
| spatialDWLS | standalone R package `spatialDWLS` (drieslab, v0.1.1, 2026-10-16; no Giotto; imports Matrix, MatrixGenerics, quadprog, stats) | Implements Dong & Yuan 2021 in two rounds: candidate cell types per cluster of spots by enrichment, then each spot re-solved with the types present. Runs on the PC, so it enters the Week 2 pilot. Inputs to pre-register: log-normalised spots; a signature built with `dwls_signature()` from the training split on the normalised (not logged) scale, with marker genes from the training split (pilot setting: top 100 genes per type by log2 fold-change of mean normalised expression, type vs the mean of the other types); **a cluster per spot** ("any clustering of the spots works"; pilot setting: k-means, k = 10, nstart 10, on the first 20 principal components of the log-normalised spots over the marker genes); `n_cell` = the true cells-per-spot level of each condition (§2.2). No coordinates needed. Limit: a new re-implementation, not the Giotto version (`runDWLSDeconv`) evaluated in Li 2022 and Spotless; version recorded in `renv.lock`. Built-in spot filter: to check in the source |
| MuSiC | MuSiC (R; v1.0.0, GitHub xuranw/MuSiC commit f21fe67, installed 2026-10-29) | Second baseline next to NNLS (DECISIONS 2026-10-26): a bulk method, among the least stable when the reference changes in Spotless. No coordinates (spots = columns of `bulk.mtx`). Needs per-cell donor labels (`samples`) from the training split. Built-in filter: genes with zero mean in the spots are dropped |
| NNLS | baseline (R package `nnls` v1.6, CRAN) | |

- **Not included:** CARD, because it exploits spatial structure that our pseudo-spots by construction do not have; with artificial coordinates it would be penalised unfairly. It is not excluded for performance: it was among the best methods in Li et al. 2023 (Nat Commun, doi:10.1038/s41467-023-37168-7). SPOTlight removed (optional in PROTOCOLLO §5.4). DECISIONS 2026-10-26.
- **Baselines in the rankings:** all five methods are ranked (DECISIONS 2026-10-30, Luca), so the baselines (NNLS, MuSiC) show whether the dedicated methods beat a simple approach (in Spotless more than half did not beat NNLS).
- **Oracle information:** the true number of cells per spot is given to every method that accepts it (spatialDWLS `n_cell`, cell2location `N_cells_per_location`); declared as oracle information, which real data do not provide (DECISIONS 2026-10-30, Luca).

- Every method receives the same reference: the scRNA training split, identical for all arms.
- Dispersion φ_g: `glmGamPoi` v1.24.0 (Bioconductor; §2.3). spatialDWLS installed 2026-10-29 = v0.1.1, GitHub drieslab/spatialDWLS commit 45d201d.
- Package versions frozen in `renv.lock` (and the Python environment for cell2location: `[LUCA]`).
- Default parameters unless stated here, as in Crowell et al. 2023 ("methods were run using default parameters and, if available, following recommendations given by the authors"); any tuning is pre-registered: `[PILOT]`.
- Built-in spot filters kept at their defaults (checked in the package source on GitHub, `master`; to re-check against the versions in `renv.lock`): RCTD drops spots with fewer than 100 UMI (`create.RCTD(UMI_min = 100)`, spacexr `R/spacexr.R`); cell2location: no spot filter (only a gene filter on the reference, `filter_genes(cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)` in the official tutorial; checked 2026-10-23); NNLS: none; MuSiC: genes with zero mean in the spots; spatialDWLS: to check. RCTD also needs at least 25 reference cells per cell type (`CELL_MIN_INSTANCE = 25`) and integer counts (`require_int = TRUE`). Spots dropped this way are failures (7.2).
- Fallback without VPN by 30/10: ADEMP approved with a CPU pilot (DECISIONS 2026-10-14 named RCTD, CARD, NNLS; with the new method list: RCTD, spatialDWLS, MuSiC, NNLS).

## 5. Performance measures

### 5.1 Aim 1 — fidelity (per arm × tissue × cell type)
- Correlation of log-means between source and held-out (PROTOCOLLO §5.5). Gene list: the §2.6 list chosen on the training split. Pearson correlation on the log-profiles, Spearman as sensitivity (DECISIONS 2026-10-30, Luca); pseudocount `[OPEN]`. Computed on each of the 4 donors of the held-out group (§2.5); E′ on all 6 pairs of the group.
- **Noise ceiling** (DECISIONS 2026-10-26): for each tissue and cell type, the cells of the held-out donor are split at random into two halves of equal size within each cell type (n_k/2 cells each; if n_k is odd, the extra cell goes to one of the two halves at random), so that Spearman–Brown's equal-halves assumption holds; the correlation r between the two half-profiles (same correlation, same genes as above) is corrected with Spearman–Brown, reliability = 2r / (1 + r); √reliability is the highest correlation any source can reach with the held-out. Repeated over B random splits (B and seed: `[PILOT]`), averaged. Cell types with few held-out cells get a lower ceiling; fidelity is read against it.
- Overlap of the top-k markers (k: `[LUCA]`).
- Cell-type separability (metric: `[LUCA]`).
- Monte Carlo SE only where μ is random (P, M, C−, possibly G); E and E′ are fixed given the split. `[LUCA]`

### 5.2 Aim 2 — error at two levels (Morris §5.2, adapted)
- **Inside replicate i:** RMSE_i = √(mean over spots × cell types of (π̂_sk − π_sk)²); JSD_i = mean over spots of the Jensen–Shannon divergence between π̂_s and π_s.
- **Across replicates:** mean RMSE over n_sim replicates, with MCSE = SD(RMSE_i) / √n_sim; same for JSD.
- Per method × arm × condition, and on the anchor.
- This per-replicate RMSE over a vector of proportions is our adaptation; Morris defines MSE for one scalar θ.

### 5.3 Aim 2 — rankings
- Rank the methods by mean RMSE (primary) and by mean JSD (secondary) within each arm × condition.
- **Kendall τ** between each arm's ranking and the anchor ranking.
- **Kendall W** across arms (agreement of all arms).
- MCSE of τ and W: bootstrap over replicates (B: `[PILOT]`).
- With 5 methods there are only 10 pairs, so τ takes 11 values (−1 to 1 in steps of 0.2); a one-swap change is τ = 0.8. With 3 ranked methods (if the baselines are left out, section 4) there are 3 pairs and one swap gives τ = 0.33. Report also the pairwise RMSE differences between methods with their MCSE, so a swap within Monte Carlo error is not read as a real change (rule in 5.4).

### 5.4 Decision rules
- A ranking **differs** from the anchor (or from another arm) if at least one difference between two methods changes sign beyond 2 × MCSE (DECISIONS 2026-10-30, Luca).
- An arm is **closest** to the anchor: proposal, highest Kendall τ with the anchor, ties broken by the mean |ΔRMSE| against the anchor `[OPEN]`.

## 6. Number of replicates (n_sim)

- Key measure: mean RMSE per method × arm × condition `[LUCA]` (or the between-method difference).
- n_sim = (SD_upper / target MCSE)², where SD_upper is the upper 95% confidence limit of SD(RMSE_i) from the pilot (chi-square, 9 df: about 1.83 × SD with 10 replicates), because with 10 replicates the SD itself is uncertain by about 24% (DECISIONS 2026-10-30, Luca) `[PILOT]`; target MCSE `[OPEN]`. Halving the MCSE needs 4× the replicates.
- Take the largest SD across the cells of the design, so every cell reaches the target.
- Check the cost: n_sim × DGMs × methods × runtime per run `[PILOT]`; cell2location may get its own n_sim by amendment.

## 7. Checking the estimates: failures (Morris §5.1)

### 7.1 No extra filter on the pseudo-spots
PROTOCOLLO §5.2 has no quality filter on the pseudo-spots, and we add none: the real data are filtered upstream (2.6, if Luca approves), and each method applies only its own built-in filter (section 4), as a user running the published pipeline would. Spots that a method drops or cannot handle are failures (7.2).
- Depth levels (2.2) are 10,000 and 20,000 UMI per spot, far above the built-in threshold of RCTD (100 UMI per spot). Still to check: any threshold of spatialDWLS. The general rule stays: check where the low level falls relative to the built-in threshold of RCTD (100 UMI per spot) and any threshold of spatialDWLS. If it is close, many spots would be dropped and the difference between methods would come from the filter rather than from the deconvolution. `[REAL DATA]`

### 7.2 Failures
1. **Count** failures per method × arm × condition (non-convergence, dropped spots, proportions not summing to 1 within tolerance `[LUCA]`) and report them as a result. Never drop them silently: failures cluster in the hardest datasets, so dropping them makes a method look better.
2. **Main analysis:** rank on the spots returned by **all** methods in that replicate.
3. **Sensitivity analysis:** failed spots get the worst score; recompute RMSE, τ and W and report whether the ranking changes.
- Failures can still differ by arm through μ (e.g. C−). Which failures actually occur: `[PILOT]`.

## 8. Reproducibility and seeds

- One seed per study `[set at approval]`; store the RNG state at the start of each replicate.
- **Check that deconvolution packages do not reset the seed internally:** compare the RNG state before and after two runs. If they do, replicates are correlated → fix before the main run.
- Parallel runs: one stream per process (L'Ecuyer-CMRG / rstream), derived from the single study seed.
- Save the simulated datasets (or their seeds and states) so any replicate can be regenerated.

## 9. Counter-proofs (PROTOCOLLO §5.6)

1. **C−** must give low fidelity and rankings different from the anchor. If rankings do not change even with meaningless profiles, the benchmark is insensitive to expression level: that is a result and is reported.
2. **E′** is the reference for donor-to-donor variability, not a ceiling: a synthetic source (P, M or G) above E′ is not evidence of leakage and is reported with its explanation (for M or G: a population-average profile, estimated from many donors, can be closer to a new donor than a single other donor is; P, built de novo, has no such explanation, so a P result above E′ is examined further). The ceiling is the noise ceiling of 5.1; a source above it beyond Monte Carlo error points to an error or leakage and is investigated before going on. *(Replaces PROTOCOLLO §5.6 point 2: DECISIONS 2026-10-26.)*
3. **Leakage, checked on the provenance of the data** (DECISIONS 2026-10-26), always and not only when a source looks too good: (a) automatic test that no donor ID appears in two splits (for colon, both samples of a donor in the same split); (b) for arm G, only held-out and E′ donors whose processed data became public after July 2024 (GEM-1 cutoff, section 2.1).

## 10. Pilot (Week 2)

Purpose: fill the `[PILOT]` placeholders, not to look at results.
- Scope (Pierangela 2026-10-29, approved by Luca 2026-10-30): **all three tissues**, run in the order colon (to debug), lymphoid, brain; arms E and C−, plus the real anchor (added 2026-10-30, because n_sim takes the largest SD across the design and the anchor is part of every condition); methods RCTD, spatialDWLS, MuSiC, NNLS (CPU); **10 replicates** per arm × condition (2 cells-per-spot × 2 depth levels). Why all three: n_sim takes the largest SD(RMSE_i) across the design (§6), and the SD depends on the biology of each tissue (19, 27 and 30 reference types; how similar the types are; how many held-out cells the anchor has); only a brain pilot can show whether the compute fallback of §2.5 is needed. With 10 replicates the relative SE of the SD is about 1/√(2·9) ≈ 24%. If brain is too slow on the PC, it gets fewer replicates and that measurement is the evidence for the fallback.
- **Donors (approved by Luca 2026-10-30):** the pilot never uses the held-out or E′ donors. In each tissue a pilot held-out and a pilot E′ donor are drawn from the training split with `assign_split()` and a separate seed (20261030); the pilot references are built from the remaining training donors. The settings fixed by the pilot (n_sim, number of spots) are therefore chosen inside training, as with a validation split in machine learning, and the test donor is used once, in the main run.
- No extra datasets: the pilot uses the three datasets of §2.5 (PROTOCOLLO §6 asks only for times, cost and variability).
- Outputs used: SD(RMSE_i) for n_sim, runtime per run, failure types, seed check.
- Pilot results do not enter the thesis results.

## 11. Declared limits (for the Discussion)
- Aim 1 uses a group of 4 randomly drawn held-out donors per tissue (distribution of donor-to-donor variability); Aim 2 uses one held-out donor (plus a second as sensitivity if time allows), so its rankings are conditional on that donor. Full internal-external cross-validation (each donor left out in turn; Steyerberg & Harrell 2016, doi:10.1016/j.jclinepi.2015.04.005; Debray et al. 2023, doi:10.1136/bmj-2022-071058; leave-one-subject-out in MuSiC, Wang et al. 2019) is not done for compute and time (DECISIONS 2026-10-29, 2026-10-30).
- The true number of cells per spot is given to the methods that accept it (oracle information, DECISIONS 2026-10-30): real data do not provide it.
- Lymphoid data are 10x 5′ (Wells 2025): this may disadvantage GEM-1 for technical reasons (DECISIONS 2026-10-26).
- Potential limit for G: the CELLxGENE version used to train the gem-1-sc served by the API is not known (DECISIONS 2026-10-30).
- Training sizes differ between tissues (brain 197 donors, colon 16, lymphoid 17), together with chemistry (nuclei 3′ v3, cells 3′ v2, cells 5′ v2): patterns compared between tissues confound tissue with training size and technology.
- spatialDWLS is run through a new standalone re-implementation, not the Giotto version evaluated in earlier benchmarks.
- Pseudo-spots have no geometry (PROTOCOLLO §8), which is why CARD is not included.

## 12. Amendments

| Date | Section | Change | Reason | Approved by |
|---|---|---|---|---|
| | | | | |
