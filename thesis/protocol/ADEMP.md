# ADEMP protocol — expression sources and deconvolution benchmarks

**Status:** DRAFT v0.1, 2026-10-26 (v0 of 2026-10-20 updated with Luca's answers, DECISIONS 2026-10-26). Not approved, not pre-registered. Becomes valid only when Luca approves it with a date (Week 2 checkpoint).
**Sources:** `PROTOCOLLO.md` §1, §2, §5; `DECISIONS.md`; Morris, White & Crowther 2019 (doi:10.1002/sim.8086), notes in `writing/letture.md`.
**Placeholders:** `[PILOT]` = fixed from the Week 2 pilot · `[REAL DATA]` = fixed from the training split of the real data · `[LUCA]` = decision for Luca.
Any change after approval is a dated amendment (section 12), written before the results it could affect.

---

## 1. Aims

- **Aim 1 — Fidelity.** How close are the cell-type mean profiles μ (genes × cell types) of each expression source to real **held-out** profiles of the same tissue and cell type?
- **Aim 2 — Benchmark sensitivity.** Do the error and the ranking of deconvolution methods change with the expression source? Which source gives the ranking closest to the one obtained on pseudo-spots built from **real held-out cells** (the real anchor)?

### Hypotheses `[LUCA]`
- **H1 (negative control).** Arm C− has the lowest fidelity in every tissue, and its method ranking differs from the anchor.
- **H2 (donor-to-donor reference).** Arm E′ shows how close a second real donor gets to the held-out donor, in fidelity and in ranking agreement with the anchor. It gives the scale on which the other arms are read: a source close to E′ is "as good as another real donor". It is not a ceiling. A synthetic arm may reach or exceed E′ without anything being wrong: M (fitted on the training donors) and G (trained on public archives) are estimated from many donors, so they give a population-average profile, which can be closer to a new donor than one other donor is. Exceeding E′ is therefore not evidence of leakage. The ceiling of Aim 1 is the split-half noise ceiling of the held-out donor (5.1); leakage is checked on the provenance of the data (section 9). *(DECISIONS 2026-10-26: E′ as donor-to-donor reference, Pierangela, approved by Luca; noise ceiling and provenance-based leakage check, Luca.)*
- **H3 (main question).** At least one synthetic arm (P, M, G) gives a method ranking that differs from the anchor ranking.
- Thresholds that turn "differs" and "closest" into decisions: `[LUCA]` (see section 5.4).

## 2. Data-generating mechanisms

### 2.1 Main factor: source of μ
| Arm | Source of μ | Role |
|---|---|---|
| E | Empirical means from the real scRNA training split | "Honest" reference |
| P | Parametric de novo: means from a gamma distribution, markers with controlled fold-changes | Maximum control, minimum realism |
| M | scDesign3 trained on the training split | Model-based middle ground |
| G | GEM-1-sc, "mean estimation" mode, from CL/UBERON metadata | Generative model (optional, gate G0 by 30/10) |
| C− | E profiles with gene labels permuted within each cell type | Negative control: must fail |
| E′ | Means from a second, independent real donor | Reference for donor-to-donor variability (not a ceiling) |

- P parameters (gamma shape/rate, number of markers, fold-change): `[REAL DATA]` (matched to the E arm's mean–variance range).
- C− permutation: one new permutation per replicate, or one fixed per study? Proposal: one per replicate, so the result does not depend on one particular shuffle. `[LUCA]`
- G is evaluated only in tissues whose held-out and E′ donors are post-cutoff and come from the same dataset as the other arms: brain (PsychAD) and lymphoid (Wells 2025). Colon (Pelka 2021) is so far "not evaluated" for G (DECISIONS 2026-10-26), but must be rechecked: Pelka appears on CELLxGENE only since 2025-08-05 (`data_inventory.xlsx`), so it may be post-cutoff under the criterion below. `[TO VERIFY]` Post-cutoff = not in the July 2024 version of CELLxGENE, the single-cell training data of GEM-1 (preprint, doi:10.1101/2025.09.08.674753, Methods; to be confirmed for the model served by the API at gate G0).
- If G0 fails, G is dropped and the design has 5 arms. Cut rule if late: arm G, then the brain tissue (DECISIONS 2026-10-26).

### 2.2 Secondary factors (fully crossed, fixed grid)
| Factor | Levels |
|---|---|
| Tissue | colon/CRC, lymphoid, brain |
| Cells per spot | 2 levels: `[REAL DATA]` |
| Sequencing depth | 2 levels: `[REAL DATA]` |

Full design: 5 arms (E, P, M, C−, E′) × 3 tissues × 2 × 2 = 60 DGMs, plus arm G on 2 tissues × 2 × 2 = 8 → **68 DGMs** (72 if colon also qualifies for G; 60 without G), plus the real anchor for each tissue × cells-per-spot × depth (12 conditions).
- Number of cell types per tissue and their labels: `[REAL DATA]`.
- Number of spots per pseudo-dataset: `[PILOT]`.
- Dirichlet concentration α (how mixed spots are): `[REAL DATA]` / `[LUCA]`.

### 2.3 Generating one replicate
1. For each spot, draw a composition from a Dirichlet(α) and sample the cells.
2. For each cell of type k, the count of gene g is negative binomial with mean μ_gk (scaled by library size) and dispersion φ_g estimated from the real training data. **The noise model is the same for all arms**; only μ changes.
3. Sum the cells of the spot.

**Real anchor:** pseudo-spots built by summing **real cells** of the held-out donor, with the same Dirichlet compositions. The truth is known because each cell's type is known.
- Colon: every pseudo-spot is built inside a single sample (normal or tumour), never mixing the two, because such a spot does not exist (DECISIONS 2026-10-26). Whether the simulated arms use one μ for normal + tumour or two separate ones depends on the Pelka reference choice: `[OPEN]`.
- How the depth factor applies to the anchor. Proposal: real cells already have their own counts, so downsample the anchor spots at random to the same target depth as the simulated arms. `[LUCA]`

### 2.4 Common random numbers
This follows from PROTOCOLLO §5.2 (only the source of μ changes; the anchor uses "la stessa composizione Dirichlet"). Within replicate i, all arms and the anchor use the **same** compositions and the same cell-type draws per spot; only μ (and therefore the counts) differs. Differences between arms are then not inflated by different compositions, which lowers the Monte Carlo SE of between-arm comparisons (Morris §5.3). All methods are run on the same datasets (Morris recommendation).

### 2.5 Data split
Split by donor, never by cell (DECISIONS 2026-10-14 and 2026-10-26). One dataset per tissue; training, held-out and E′ are different donors of the same dataset (DECISIONS 2026-10-26). Training split → E, M, φ_g, P calibration and the reference given to every method. Held-out donor → Aim 1 target, noise ceiling and the anchor. Second independent donor → E′.
- Colon: Pelka donors have a normal and a tumour sample; both samples of a donor always go to the same split (no pseudoreplication, no leakage).

| Tissue | Dataset | Chemistry | Donors: training | held-out | E′ |
|---|---|---|---|---|---|
| Brain | PsychAD (DLPFC) | 10x 3′ v3, single-nucleus | `[REAL DATA]` | `[REAL DATA]` | `[REAL DATA]` |
| Colon/CRC | Pelka 2021, 3′ v2 subset | 10x 3′ v2 | `[REAL DATA]` | `[REAL DATA]` | `[REAL DATA]` |
| Lymphoid | Wells 2025, lymph node | 10x 5′ v2 | `[REAL DATA]` | `[REAL DATA]` | `[REAL DATA]` |

Donor counts per split are filled by `scripts/02_split.R` (Luca's request, 2026-10-26); candidates and rules in `data_inventory.xlsx`, sheet `selected`.

### 2.6 Filtering of the real data (before building μ) — proposal `[LUCA]`
PROTOCOLLO §5.2 does not mention this step; it is proposed to Luca. The same minimal filter is applied to every real dataset before μ, φ_g, the reference for the methods and the anchor are computed. It follows Crowell et al. 2023 (Additional file 1, §6.1.2), who filtered their real references "as insufficient numbers of cells and low-quality observations can interfere with estimation of simulation parameters". Filtered data are saved in `data/derived/`; the downloaded files are never changed (PROTOCOLLO §3.7).
1. **Cell types:** remove a cell type from a split (training, held-out donor, E′ donor) if it has fewer than 50 cells there. In Crowell the groups are batches or clusters; for us they are cell types. The cell types removed are listed in the report.
2. **Genes:** keep genes with a count greater than 1 in at least 10 cells. Proposal: decide the gene list on the training split only and use the same list for every split and arm, so that all μ matrices have the same genes and the held-out data do not influence the choice.
3. **Cells:** keep cells with at least 100 detected genes (non-zero counts).

## 3. Estimands

- **Primary (Aim 2):** the proportion π_sk of each cell type k in each spot s. It is a vector per spot, not a single scalar θ as in Morris; we score it with RMSE and Jensen–Shannon divergence (section 5).
- **Ranking target (Aim 2):** the ordering of the methods by mean error within each arm × condition, compared with the ordering on the anchor.
- **Secondary (Aim 1):** the cell-type profiles μ and their markers, compared with the held-out real profiles.

## 4. Methods

| Method | Package / language | Notes |
|---|---|---|
| RCTD | spacexr (R) | mode: `[LUCA]`. Cable et al. 2022 recommend doublet mode only when pixels rarely hold three or more cell types (Slide-seq), otherwise full mode (e.g. Visium); our pseudo-spots mix several cell types, so the proposal is `doublet_mode = "full"` (package default is "doublet", so it must be set). No coordinates needed: `SpatialRNA(..., use_fake_coords = TRUE)` (checked in the source, 2026-10-23). See `letture.md` |
| cell2location | Python, GPU on the server | n_sim may be lower; fixed by dated amendment in the first days on the server (DECISIONS 2026-10-14). Needs raw counts. Hyper-parameters: `N_cells_per_location` set to the true cells-per-spot level of each condition (§2.2) `[LUCA]`; `detection_alpha` 20 or 200 `[LUCA]` |
| spatialDWLS | standalone R package `spatialDWLS` (drieslab, v0.1.1, 2026-10-16; no Giotto; imports Matrix, MatrixGenerics, quadprog, stats) | Implements Dong & Yuan 2021 in two rounds: candidate cell types per cluster of spots by enrichment, then each spot re-solved with the types present. Runs on the PC, so it enters the Week 2 pilot. Inputs to pre-register: log-normalised spots; a signature built with `dwls_signature()` from the training split on the normalised (not logged) scale, with marker genes from the training split `[PILOT]`; **a cluster per spot** ("any clustering of the spots works", e.g. Leiden on expression: method and resolution `[PILOT]`); `n_cell` = the true cells-per-spot level of each condition (§2.2). No coordinates needed. Limit: a new re-implementation, not the Giotto version (`runDWLSDeconv`) evaluated in Li 2022 and Spotless; version recorded in `renv.lock`. Built-in spot filter: to check in the source |
| MuSiC | MuSiC (R) | Second baseline next to NNLS (DECISIONS 2026-10-26): a bulk method, among the least stable when the reference changes in Spotless. No coordinates (spots = columns of `bulk.mtx`). Needs per-cell donor labels (`samples`) from the training split. Built-in filter: genes with zero mean in the spots are dropped |
| NNLS | baseline (R) | |

- **Not included:** CARD, because it exploits spatial structure that our pseudo-spots by construction do not have; with artificial coordinates it would be penalised unfairly. It is not excluded for performance: it was among the best methods in Li et al. 2023 (Nat Commun, doi:10.1038/s41467-023-37168-7). SPOTlight removed (optional in PROTOCOLLO §5.4). DECISIONS 2026-10-26.
- **Baselines in the rankings:** proposal: rank all five methods, so the baselines (NNLS, MuSiC) show whether the dedicated methods beat a simple approach (in Spotless more than half did not beat NNLS). Alternative: rank only RCTD, cell2location and spatialDWLS and report the baselines alongside. `[OPEN]`

- Every method receives the same reference: the scRNA training split, identical for all arms.
- Package versions frozen in `renv.lock` (and the Python environment for cell2location: `[LUCA]`).
- Default parameters unless stated here, as in Crowell et al. 2023 ("methods were run using default parameters and, if available, following recommendations given by the authors"); any tuning is pre-registered: `[PILOT]`.
- Built-in spot filters kept at their defaults (checked in the package source on GitHub, `master`; to re-check against the versions in `renv.lock`): RCTD drops spots with fewer than 100 UMI (`create.RCTD(UMI_min = 100)`, spacexr `R/spacexr.R`); cell2location: no spot filter (only a gene filter on the reference, `filter_genes(cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)` in the official tutorial; checked 2026-10-23); NNLS: none; MuSiC: genes with zero mean in the spots; spatialDWLS: to check. RCTD also needs at least 25 reference cells per cell type (`CELL_MIN_INSTANCE = 25`) and integer counts (`require_int = TRUE`). Spots dropped this way are failures (7.2).
- Fallback without VPN by 30/10: ADEMP approved with a CPU pilot (DECISIONS 2026-10-14 named RCTD, CARD, NNLS; with the new method list: RCTD, spatialDWLS, MuSiC, NNLS).

## 5. Performance measures

### 5.1 Aim 1 — fidelity (per arm × tissue × cell type)
- Correlation of log-means between source and held-out (PROTOCOLLO §5.5). Gene list: the §2.6 list chosen on the training split; Pearson or Spearman and the pseudocount: `[OPEN]`.
- **Noise ceiling** (DECISIONS 2026-10-26): for each tissue and cell type, the cells of the held-out donor are split at random into two halves; the correlation r between the two half-profiles (same correlation, same genes as above) is corrected with Spearman–Brown, reliability = 2r / (1 + r); √reliability is the highest correlation any source can reach with the held-out. Repeated over B random splits (B and seed: `[PILOT]`), averaged. Cell types with few held-out cells get a lower ceiling; fidelity is read against it.
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
- Caveat: with 5 methods there are only 10 pairs, so τ takes 11 values (−1 to 1 in steps of 0.2); a one-swap change is τ = 0.8. With 3 ranked methods (if the baselines are left out, section 4) there are 3 pairs and one swap gives τ = 0.33. Report also the pairwise RMSE differences between methods with their MCSE, so a swap within Monte Carlo error is not read as a real change. `[LUCA]`

### 5.4 Decision rules `[LUCA]`
- When does a ranking "differ" from the anchor? (e.g. τ below a threshold, or a pairwise difference that changes sign beyond 2 × MCSE)
- When is an arm "closest" to the anchor? (e.g. highest τ, ties broken by mean |ΔRMSE| vs anchor)

## 6. Number of replicates (n_sim)

- Key measure: mean RMSE per method × arm × condition `[LUCA]` (or the between-method difference).
- n_sim = (SD(RMSE_i) / target MCSE)², with SD from the pilot `[PILOT]`, target MCSE `[LUCA]`. Halving the MCSE needs 4× the replicates.
- Take the largest SD across the cells of the design, so every cell reaches the target.
- Check the cost: n_sim × DGMs × methods × runtime per run `[PILOT]`; cell2location may get its own n_sim by amendment.

## 7. Checking the estimates: failures (Morris §5.1)

### 7.1 No extra filter on the pseudo-spots
PROTOCOLLO §5.2 has no quality filter on the pseudo-spots, and we add none: the real data are filtered upstream (2.6, if Luca approves), and each method applies only its own built-in filter (section 4), as a user running the published pipeline would. Spots that a method drops or cannot handle are failures (7.2).
- When the two depth levels are fixed from the real data (2.2), check where the low level falls relative to the built-in threshold of RCTD (100 UMI per spot) and any threshold of spatialDWLS. If it is close, many spots would be dropped and the difference between methods would come from the filter rather than from the deconvolution. `[REAL DATA]`

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
3. **Leakage, checked on the provenance of the data** (DECISIONS 2026-10-26), always and not only when a source looks too good: (a) automatic test that no donor ID appears in two splits (for colon, both samples of a donor in the same split); (b) for arm G, only held-out and E′ donors not in the July 2024 version of CELLxGENE (GEM-1 cutoff, section 2.1).

## 10. Pilot (Week 2)

Purpose: fill the `[PILOT]` placeholders, not to look at results.
- Scope: one tissue `[LUCA]`, arms E and C−, methods RCTD, spatialDWLS, MuSiC, NNLS (CPU), small n `[LUCA]`.
- Outputs used: SD(RMSE_i) for n_sim, runtime per run, failure types, seed check.
- Pilot results do not enter the thesis results.

## 11. Declared limits (for the Discussion)
- Lymphoid data are 10x 5′ (Wells 2025): this may disadvantage GEM-1 for technical reasons (DECISIONS 2026-10-26).
- spatialDWLS is run through a new standalone re-implementation, not the Giotto version evaluated in earlier benchmarks.
- Pseudo-spots have no geometry (PROTOCOLLO §8), which is why CARD is not included.

## 12. Amendments

| Date | Section | Change | Reason | Approved by |
|---|---|---|---|---|
| | | | | |
