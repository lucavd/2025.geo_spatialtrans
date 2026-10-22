# ADEMP protocol — expression sources and deconvolution benchmarks

**Status:** DRAFT v0, 2026-10-20. Not approved, not pre-registered. Becomes valid only when Luca approves it with a date (Week 2 checkpoint).
**Sources:** `PROTOCOLLO.md` §1, §2, §5; `DECISIONS.md`; Morris, White & Crowther 2019 (doi:10.1002/sim.8086), notes in `writing/letture.md`.
**Placeholders:** `[PILOT]` = fixed from the Week 2 pilot · `[REAL DATA]` = fixed from the training split of the real data · `[LUCA]` = decision for Luca.
Any change after approval is a dated amendment (section 9), written before the results it could affect.

---

## 1. Aims

- **Aim 1 — Fidelity.** How close are the cell-type mean profiles μ (genes × cell types) of each expression source to real **held-out** profiles of the same tissue and cell type?
- **Aim 2 — Benchmark sensitivity.** Do the error and the ranking of deconvolution methods change with the expression source? Which source gives the ranking closest to the one obtained on pseudo-spots built from **real held-out cells** (the real anchor)?

### Hypotheses `[LUCA]`
- **H1 (negative control).** Arm C− has the lowest fidelity in every tissue, and its method ranking differs from the anchor.
- **H2 (positive control).** Arm E′ has the highest fidelity and the ranking closest to the anchor; no synthetic arm exceeds it beyond Monte Carlo error.
- **H3 (main question).** At least one synthetic arm (P, M, G) gives a method ranking that differs from the anchor ranking.
- Thresholds that turn "differs" and "closest" into decisions: `[LUCA]` (see section 6.4).

## 2. Data-generating mechanisms

### 2.1 Main factor: source of μ
| Arm | Source of μ | Role |
|---|---|---|
| E | Empirical means from the real scRNA training split | "Honest" reference |
| P | Parametric de novo: means from a gamma distribution, markers with controlled fold-changes | Maximum control, minimum realism |
| M | scDesign3 trained on the training split | Model-based middle ground |
| G | GEM-1-sc, "mean estimation" mode, from CL/UBERON metadata | Generative model (optional, gate G0 by 30/10) |
| C− | E profiles with gene labels permuted within each cell type | Negative control: must fail |
| E′ | Means from a second, independent real donor | Positive control: the attainable ceiling |

- P parameters (gamma shape/rate, number of markers, fold-change): `[REAL DATA]` (matched to the E arm's mean–variance range).
- C− permutation: one new permutation per replicate, or one fixed per study? Proposal: one per replicate, so the result does not depend on one particular shuffle. `[LUCA]`
- If G0 fails, G is dropped and the design has 5 arms (cut rule, DECISIONS 2026-10-14).

### 2.2 Secondary factors (fully crossed, fixed grid)
| Factor | Levels |
|---|---|
| Tissue | colon/CRC, lymphoid, brain |
| Cells per spot | 2 levels: `[REAL DATA]` |
| Sequencing depth | 2 levels: `[REAL DATA]` |

Full design: 6 arms × 3 tissues × 2 × 2 = 72 DGMs (60 without G), plus the real anchor for each tissue × cells-per-spot × depth (12 conditions).
- Number of cell types per tissue and their labels: `[REAL DATA]`.
- Number of spots per pseudo-dataset: `[PILOT]`.
- Dirichlet concentration α (how mixed spots are): `[REAL DATA]` / `[LUCA]`.

### 2.3 Generating one replicate
1. For each spot, draw a composition from a Dirichlet(α) and sample the cells.
2. For each cell of type k, the count of gene g is negative binomial with mean μ_gk (scaled by library size) and dispersion φ_g estimated from the real training data. **The noise model is the same for all arms**; only μ changes.
3. Sum the cells of the spot.

**Real anchor:** pseudo-spots built by summing **real cells** of the held-out donor, with the same Dirichlet compositions. The truth is known because each cell's type is known.
- How the depth factor applies to the anchor. Proposal: real cells already have their own counts, so downsample the anchor spots at random to the same target depth as the simulated arms. `[LUCA]`

### 2.4 Common random numbers
This follows from PROTOCOLLO §5.2 (only the source of μ changes; the anchor uses "la stessa composizione Dirichlet"). Within replicate i, all arms and the anchor use the **same** compositions and the same cell-type draws per spot; only μ (and therefore the counts) differs. Differences between arms are then not inflated by different compositions, which lowers the Monte Carlo SE of between-arm comparisons (Morris §5.3). All methods are run on the same datasets (Morris recommendation).

### 2.5 Data split
Split by donor or sample, never by cell (DECISIONS 2026-10-14). Training split → E, M, φ_g, P calibration and the reference given to every method. Held-out donor → Aim 1 target and the anchor. Second independent donor → E′. Donor assignment per tissue: `[REAL DATA]` (see `data_inventory.xlsx`, sheet `selected`).

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
| RCTD | spacexr (R) | mode: `[LUCA]` (doublet vs full) |
| CARD | CARD (R) | uses spatial coordinates: pseudo-spots have none → setting `[LUCA]` |
| cell2location | Python, GPU on the server | n_sim may be lower; fixed by dated amendment in the first days on the server (DECISIONS 2026-10-14) |
| NNLS | baseline (R) | |
| SPOTlight | R, optional | first to be cut if late |

- Every method receives the same reference: the scRNA training split, identical for all arms.
- Package versions frozen in `renv.lock` (and the Python environment for cell2location: `[LUCA]`).
- Default parameters unless stated here, as in Crowell et al. 2023 ("methods were run using default parameters and, if available, following recommendations given by the authors"); any tuning is pre-registered: `[PILOT]`.
- Built-in spot filters kept at their defaults (checked in the package source on GitHub, `master`; to re-check against the versions in `renv.lock`): RCTD drops spots with fewer than 100 UMI (`create.RCTD(UMI_min = 100)`, spacexr `R/spacexr.R`); CARD drops spots with fewer than 100 total counts and genes detected in 5 spots or fewer (`createCARDObject`, CARD `R/utilities.R`); cell2location: to check; NNLS: none. Spots dropped this way are failures (7.2).
- Fallback without VPN by 30/10: ADEMP approved with a CPU pilot (RCTD, CARD, NNLS).

## 5. Performance measures

### 5.1 Aim 1 — fidelity (per arm × tissue × cell type)
- Correlation of log-means between source and held-out.
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
- Caveat: with 4 methods there are only 6 pairs, so τ takes 7 values (−1 to 1 in steps of 1/3). A one-swap change is τ = 0.67. Report also the pairwise RMSE differences between methods with their MCSE, so a swap within Monte Carlo error is not read as a real change. `[LUCA]`

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
- When the two depth levels are fixed from the real data (2.2), check where the low level falls relative to the built-in thresholds of RCTD and CARD (100 UMI per spot). If it is close, many spots would be dropped and the difference between methods would come from the filter rather than from the deconvolution. `[REAL DATA]`

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
2. **E′** is the ceiling. A synthetic source that exceeds it suggests leakage (same donor in training and test, or GEM-1 training data overlapping ours) and is investigated before going on.
3. **Donor split** checked by an automatic test: no donor ID appears in two splits.

## 10. Pilot (Week 2)

Purpose: fill the `[PILOT]` placeholders, not to look at results.
- Scope: one tissue `[LUCA]`, arms E and C−, methods RCTD, CARD, NNLS (CPU), small n `[LUCA]`.
- Outputs used: SD(RMSE_i) for n_sim, runtime per run, failure types, seed check.
- Pilot results do not enter the thesis results.

## 11. Amendments

| Date | Section | Change | Reason | Approved by |
|---|---|---|---|---|
| | | | | |
