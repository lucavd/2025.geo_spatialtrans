# Readings — five lines each (PROTOCOLLO.md §10)

Format: question · data · main result · limit · what we need from it.

---

## Koytiger, Walsh et al. 2025 — GEM-1 preprint
*Generative genomics accurately predicts future experimental results.* bioRxiv 2025. doi:10.1101/2025.09.08.674753. Not peer reviewed; every author is a Synthesize Bio employee or founder.

1. **Question:** Can a generative model trained on public human transcriptomes produce realistic expression profiles from metadata alone (tissue, cell type, perturbation, technology)?
2. **Data:** 
    - **(a) TRAINING**
    - Bulk: human Illumina RNA-seq in SRA up to 30 Jun 2024 (470,691 samples after QC)
    - Single-cell: CZI CELLxGENE (including Tabula sapiens v1) up to July 2024 plus 10 scPerturb datasets, about 41.5 M cells, 669 cell types, 259 tissues
    Metadata harmonised by an LLM agent and mapped to ontologies (Cell Ontology, Uberon, …)
    - **(b) TEST**
    - Bulk: SRA data 1 Jul – 30 Sep 2024
    - Single-cell: Tabula Sapiens v2 (=9 new donors)
3. **Main result:** In metadata-only mode, accuracy reaches the level of pseudoreplicates (same metadata, different labs), except for novel cell lines. Single-cell generation accuracy is fairly flat across tissues, and its embeddings annotate cell types comparably to scGPT and SCimilarity.
4. **Limit:** Closed weights and self-run evaluation by the company, so it cannot be independently checked. About 95% of the "future" test set comes from already-seen contexts, so the test is mostly interpolation. There is no simple baseline (e.g. "no change" or "average profile"). The headline metric `r_sample` stays near 0.9 even for unrelated samples. Single-cell pseudoreplicates come from a single consortium, which the authors themselves call optimistic. Single-cell evaluation uses only `r_sample` and MSE.
5. **What we need from it:** It describes arm G's source (single-cell model, mean profiles from CL/UBERON metadata). It gives a stated training cutoff: CELLxGENE July 2024 for the single-cell model, so arm G's held-out donors should be deposited after that date. This must be confirmed with the company, because the deployed model may have been retrained since (gate G0). Their "replicate/pseudoreplicate ceiling" is the same logic as our E′ arm. Their weak points (inflated whole-profile metrics, no baseline) are reasons to use stringent per-cell-type measures and the C− control in our design.

## Morris, White, Crowther 2019
*Using simulation studies to evaluate statistical methods.* Statistics in Medicine 2019;38(11):2074–2102. doi:10.1002/sim.8086. Peer reviewed. The authors declare that they deliver a short course on simulation studies, from which this paper grew and from which they have benefited financially.

1. **Question:** How should simulation studies be planned, run, analysed and reported? The paper gives a structured approach, the ADEMP model: Aims, Data-generating mechanisms, Estimands, Methods, Performance measures.
2. **Data:** Review of 100 simulation studies, randomly sampled from those in *Statistics in Medicine* vol. 34 (2015), and by a simple worked example (two data-generating mechanisms, three methods).
3. **Main result:** The paper defines the performance measures (bias, empirical SE, MSE, model SE, coverage, power), each with its Monte Carlo SE. It shows how to choose n_sim from the Monte Carlo SE you want for the key measure, and how to present results (e.g. graphically).
4. **Limit:** Some Monte Carlo SE formulas are approximate or assume normally distributed estimates. The guidance is built around a single scalar estimand (one θ, e.g. a treatment effect) and its confidence interval. For our setting, we need to adapt this to:
    - inside each benchmark (deconvolution methods RCTD, CARD, ...): the estimand is the vector of cell-type proportions of each spot, so we score it to RMSE and Jensen–Shannon divergence instead of the bias of a single θ;
    - across sources of expression (arms E, P, M, G, ...): the target is the **ranking** of the methods, compared with Kendall's τ (agreement with the real anchor) and Kendall's W (agreement across arms), each with its Monte Carlo SE.
5. **What we need from it** (for `protocol/ADEMP.md`):
    - **ADEMP skeleton:** 
        - **Aim** = do benchmark rankings change with the expression source; DGM = Dirichlet pseudo-spots with NB counts, one per arm; 
        - **Estimand** = true cell-type proportions per spot;
        - **Methods** = RCTD, CARD, cell2location, NNLS; 
        - **Performance** = RMSE and Jensen–Shannon per method × arm, plus the rankings.
    - **Performance at two levels (§5.2):** inside each replicate, one RMSE over all spots × cell types (how wrong the method was on that dataset); across replicates, the mean RMSE with its Monte Carlo SE = SD(RMSE_i) / √n_sim. Morris gives the MSE for a single scalar θ; the per-replicate RMSE over a vector of proportions is our adaptation.
    - **Number of replicates:** n_sim = (SD / target MCSE)², with the SD taken from the pilot. Halving the MCSE needs 4× the replicates.
    - **Checking the estimates data (§5.1):** count failures (non-convergence, dropped spots, proportions not summing to 1) per method × arm and report them as a result. Never drop them silently: failures cluster in the hardest datasets, so dropping them makes the method look better.
    - **Failures in our design:** the protocol (§5.2) has no quality filter on the pseudo-spots, so every spot goes to every method. The rules above apply to the spots a method fails on: count them per method × arm, rank on the spots that all methods returned, and as a check give the failed spots the worst score. Failures can differ between arms only through μ (e.g. C−), because noise and depth are the same for all arms.
    - **Seeds:** set the seed once per study and store the random-number state at the start of each repetition. **Check that the deconvolution packages do not reset the seed internally** (compare the state before and after two runs), otherwise replicates are correlated. In parallel runs, give each process its own stream (L'Ecuyer-CMRG / rstream) derived from the single seed.



## Crowell, Morillo Leonardo, Soneson, Robinson 2023 (+ Author Correction 2024)
*The shaky foundations of simulating single-cell RNA sequencing data.* Genome Biology 2023;24:62. doi:10.1186/s13059-023-02904-1. Peer reviewed; no competing interests declared. Author Correction (Genome Biol 2024;25:178, doi:10.1186/s13059-024-03329-0): only the legend of Table 1 was wrong (it described a column that does not exist); the table itself and all results are unchanged. Our PDF is the corrected version.

1. **Question:** When we simulate single-cell data, how similar are the simulated data to the real data the simulator learned from? And if we compare analysis methods (a benchmark) on simulated data, do we reach the same conclusions as on real data?
2. **Data:** Of the 65 tools tagged "simulation" in the scRNA-tools database, they kept the 16 that could be installed and that learn their parameters from a real dataset ("reference-based"). Each simulator was fitted to subsets of 12 public scRNA-seq datasets of three kinds: one group of cells (10 references), several batches (8), several clusters (8). To measure similarity they computed simple summaries on real and simulated data (e.g. mean and variance of each gene, fraction of zeros, library size, correlation between cells) and compared their distributions with the Kolmogorov–Smirnov distance. Then they ran two benchmarks on both the real and the simulated version of the same data: 6 batch-correction methods and 9 clustering methods.
3. **Main result:**
    - Most simulators can generate only one group of cells. To get several groups, the user has to add artificial differences by hand. The differences *between* groups (batches, clusters) were the worst-reproduced feature.
    - Most realistic overall: ZINB-WaVE, scDesign2 and muscat.
    - Batch correction: for most simulators the results on simulated data were fairly close to those on real data, only somewhat too optimistic.
    - Clustering: with every simulator, the methods scored better on simulated data than on real data. The ranking of the clustering methods agreed poorly with the ranking on real data. Among simulators of similar quality, the ranking depended more on which real dataset was used than on which simulator.
    - The most realistic simulator is not always the best one for a benchmark: muscat was among the most realistic, but its rankings disagreed with real data; SPsimSeq was only average in realism, but gave the batch-correction ranking closest to real data.
    - Their own summary: "The most truthful model for real data is real data."
4. **Limit:** Single-cell data only: no spatial data and no deconvolution. Only two tasks (batch correction, clustering). Realism is a plain average over all summaries, and the authors themselves say that which summaries matter depends on the task. In the main text, the agreement between rankings is a Spearman correlation without any uncertainty (I have not read their Additional file 1, which may say more). Simulators that do not learn from real data ("de novo") were left out by design.
5. **What we need from it:**
    - **Motivation:** protocol §1 cites this paper: the choice of simulator changes the conclusions of single-cell benchmarks, and nobody has measured this for spatial deconvolution. That is our question.
    - **Why two objectives:** realism and a correct ranking can disagree (muscat vs SPsimSeq). So Objective 1 (fidelity) and Objective 2 (ranking) are measured separately, as in protocol §1 and §5.5.
    - **Real anchor:** they compared results on simulated data with results on real data. We do the same with pseudo-spots made of real held-out cells (protocol §5.2), compared with Kendall's τ (§5.5).
    - **Errors, not only rankings:** Objective 2 (protocol §1) asks whether both the ranking *and the error* of the methods change with the source. Crowell found methods looking better on simulated data, so we look at whether RMSE and JSD are lower on the synthetic arms than on the real anchor.
    - **Arm P:** P is built from scratch (protocol §5.2: "massimo controllo, minimo realismo"), the kind of simulator Crowell left out. They warn that conclusions from such simulations "should be treated with caution". Our study measures how much this matters for the ranking.
    - **Arm M:** we use scDesign3. Crowell tested the earlier scDesign2 (among the most realistic), not scDesign3.
    - **Tissues:** rankings depended on the reference dataset, which agrees with the protocol's choice of three tissues and with reporting results for each tissue.
    - **Filtering:** they filtered the real references only, not the simulated data (Additional file 1, §6.1.2: groups with fewer than 50 cells removed; genes kept if their count is greater than 1 in at least 10 cells; cells kept with at least 100 detected genes), and ran the methods with their default parameters. For us: no extra filter on the pseudo-spots and each deconvolution method with its defaults, built-in filters included (ADEMP §7.1), as the protocol already implies. Using the same three rules on our real data before building μ (cell types as groups) is a proposal for Luca (ADEMP §2.6).
    - **Uncertainty:** they compared rankings without uncertainty. Our protocol (§5.5, "Incertezza") requires a Monte Carlo SE for every measure (Morris 2019), so we avoid this weakness by design.

## Cable, Murray, Zou, Goeva, Macosko, Chen, Irizarry 2022 — RCTD
*Robust decomposition of cell type mixtures in spatial transcriptomics.* Nature Biotechnology 2022;40:517–526. doi:10.1038/s41587-021-00830-w. Peer reviewed. R package `spacexr`.

1. **Question:** Can the cell-type proportions of spatial measurements that mix several cells be estimated from a single-cell reference?
2. **Data:** Mouse brain single-cell and single-nucleus references; Slide-seq and Visium data; simulated mixtures built from single cells of the reference (pairs of cells, and three or four cell types per pixel with 10,000 UMIs "to model the conditions of Visium").
3. **Main result:** The counts of each pixel are modelled as Poisson, with a mean equal to the pixel's total UMIs times a weighted mix of the reference cell-type profiles (average expression per cell type); the weights are the proportions. It recovered the true proportions of the simulated mixtures. "Doublet mode" allows at most two cell types per pixel; "full mode" allows any number, as in Visium.
4. **Limit:** Validated mainly in mouse brain, on simulations built by the authors; cell types missing from the reference are a problem.
5. **What we need from it:** It fits our pseudo-spots, which are the same kind of mixture as its simulations. Settings: `doublet_mode = "full"` (the package default is "doublet"); no coordinates needed (`SpatialRNA(..., use_fake_coords = TRUE)`); raw integer counts; at least 25 reference cells per cell type; spots and reference cells under 100 UMI are dropped (checked in the `spacexr` source, 2026-10-23).

## Ma, Zhou 2022 — CARD (no longer in the study: it exploits spatial structure that our pseudo-spots do not have; DECISIONS 2026-10-26)
*Spatially informed cell-type deconvolution for spatial transcriptomics.* Nature Biotechnology 2022;40:1349–1359. doi:10.1038/s41587-022-01273-7. Peer reviewed. R package `CARD`.

1. **Question:** Does using the fact that neighbouring spots have similar cell-type compositions improve deconvolution, also when the single-cell reference does not match the tissue well?
2. **Data:** Simulations built from real data, with a parameter that sets how spatially correlated the compositions are, five replicates per setting and six competing methods (among them RCTD, SPOTlight, cell2location, stereoscope, spatialDWLS, MuSiC); four real datasets, including mouse olfactory bulb and human pancreatic cancer.
3. **Main result:** CARD is a non-negative matrix factorization of the spot data (normalized expression) into a reference matrix B, built from the single-cell data, times the proportions V. The proportions follow a **conditional autoregressive (CAR)** model: each spot's composition is pulled towards the weighted composition of the other spots, with weights from a Gaussian kernel of the distance between spots and a spatial autocorrelation parameter that sets the strength of the pull. CARD was the most accurate method in their simulations (median RMSE 0.079), also with mismatched references, but "the advantage of CARD over the other methods shows a clear dependency on spatial correlation". The reference matrix follows MuSiC: genes present in both datasets, cell-type-informative genes (mean in a cell type at least 1.25 log-fold above the other cell types), and removal of the 1% of genes with the highest dispersion.
4. **Limit:** CARD models normalized data, not raw counts (the authors list a count model as future work). Its gain comes from spatial correlation, so with no spatial structure its advantage should disappear. The benchmark was run by the method's authors.
5. **What we need from it:**
    - **CARD without coordinates (ADEMP §4, Open question 5a):** the package **cannot run without coordinates**: `createCARDObject()` stops with "Please provide the matched spatial location data frame" when `spatial_location` is missing (checked in the source, `R/utilities.R`, 2026-10-23). Our pseudo-spots have no positions, so we must give it coordinates. How (for example a grid with spots placed at random, so that neighbours share nothing) is a design decision to discuss with Luca; the paper suggests that CARD would then behave as a non-spatial method.
    - **Built-in filters** (checked in the source, 2026-10-23; ADEMP §4 confirmed): spots are kept only with total counts between 100 and 1,000,000 (`minCountGene = 100`, applied to spots despite its name); genes are kept only if detected in more than 5 spots (`minCountSpot = 5`, applied to genes despite its name). The package documentation describes the two parameters the other way round.
    - **Replicates:** CARD's own simulations used five replicates per setting, an example of reporting variability that our Monte Carlo SE extends.

## Kleshchevnikov, Shmatko, Dann et al. 2022 — cell2location
*Cell2location maps fine-grained cell types in spatial transcriptomics.* Nature Biotechnology 2022;40:661–671. doi:10.1038/s41587-021-01139-4. Peer reviewed. Python package `cell2location` (scvi-tools, GPU).

1. **Question:** Can a Bayesian model map fine-grained cell types, including rare ones, in spatial data more sensitively than existing tools?
2. **Data:** A synthetic dataset of 2,500 locations built by combining cells of 49 mouse-brain cell types (from single-nucleus data) with different abundance patterns; paired Visium and single-nucleus data of the mouse brain; human lymph node and human gut. Compared with RCTD, SPOTlight, Seurat, stereoscope and **NNLS**.
3. **Main result:** Two steps. First, a signature for each cell type is estimated from the single-cell or single-nucleus reference, by default with a negative binomial regression that can remove batch effects. Second, the spatial counts are modelled as negative binomial, with expected counts built from the signatures times the abundance of each cell type, plus gene-specific scaling, per-location detection sensitivity and an additive background term. It estimates absolute cell abundance per location, not only proportions. It was more accurate than the other methods on the synthetic data and resolved fine immune populations in lymph node and gut.
4. **Limit:** Two hyper-parameters must be set by the user: the expected number of cells per location N̂ ("an appropriate setting of this prior is essential to inform the estimation of absolute cell type abundance values", although the model "is robust to a range of similar values"), and the regularization of detection sensitivity. It needs a GPU for large data. The synthetic benchmark was built by the authors.
5. **What we need from it:**
    - **Raw counts:** "Cell2location needs untransformed unnormalised spatial mRNA counts as input" (official tutorial).
    - **Built-in filters** (ADEMP §4 said "to check"; checked 2026-10-23): **no spot-level filter.** The tutorial removes no spots; it only removes mitochondrial genes from the spatial data and filters genes of the reference with `filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)` (the function's own defaults are 15, 0.05 and 1.12; source `cell2location/utils/filtering.py`).
    - **N̂ (`N_cells_per_location`):** the tutorial uses 30 and says to adapt it to each tissue. In our design the true number of cells per spot is known (the "cells per spot" factor, ADEMP §2.2), so N̂ can be set to the level of each condition; to fix in the ADEMP.
    - **`detection_alpha`:** 20 (tutorial default for human Visium with technical variability) or 200 (low technical variability). Our pseudo-spots have no per-spot sensitivity variation beyond library size, so 200 may fit; to decide in the ADEMP.
    - **NNLS:** cell2location's own benchmark already included NNLS as a comparison, which supports our baseline.

## Elosua-Bayes, Nieto, Mereu, Gut, Heyn 2021 — SPOTlight (no longer in the study, DECISIONS 2026-10-26)
*SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes.* Nucleic Acids Research 2021;49(9):e50. doi:10.1093/nar/gkab043. Peer reviewed. R/Bioconductor package `SPOTlight`. Optional method in the protocol (§5.4).

1. **Question:** Can a non-negative matrix factorization (NMF), started from cell-type marker genes, deconvolute spots accurately, also with small or shallow single-cell references?
2. **Data:** Synthetic mixtures of 2–8 cells taken from single-cell datasets (if a mixture had more than 25,000 UMIs it was down-sampled to 20,000 "to better simulate biological capture locations"); PBMC references from 13 single-cell technologies and at different sequencing depths; mouse brain and human pancreatic cancer spatial data. Compared with MuSiC, CIBERSORTx, DeconRNAseq, SCDC, RCTD and NMFreg, all "run with default settings".
3. **Main result:** Two steps. First, NMF is trained on the single-cell reference, with each topic seeded with the marker genes of one cell type, giving a topic profile per cell type. Second, non-negative least squares (NNLS) maps each spot onto these topics, and a second NNLS turns the topic weights into cell-type proportions; weights below a minimum contribution are set to zero. Accuracy was high also with low cell numbers or shallow references. Performance was measured as presence/absence (sensitivity, specificity) and as proportions with the **Jensen–Shannon divergence**, the same measure as our §5.5.
4. **Limit:** The results depend on user choices: the marker genes, the number of cells per cell type used for training, and the minimum contribution threshold. The benchmark was run by the method's authors.
5. **What we need from it:**
    - **No coordinates needed:** `SPOTlight()` has no coordinates argument; the spots are only the columns of the count matrix (checked in the package source, `R/SPOTlight.R`, `devel` branch, 2026-10-23). So SPOTlight can run on our pseudo-spots as designed. This is the alternative to CARD in Luca question 5a.
    - **Built-in filters:** no spot filter in the function. Genes come from the marker list (`mgs`, required) and optional highly variable genes (`hvg`). Proportions under 1% are set to zero (`min_prop = 0.01`).
    - **Extra choice it brings:** the marker genes must be computed by us on the training split (the paper used Seurat's `FindAllMarkers`); this must be fixed in the ADEMP if SPOTlight is used.
    - **Our generator:** their synthetic mixtures (a few cells summed, then down-sampled to a fixed depth) are close to our pseudo-spots with the depth factor.

## Song, Wang, Yan, Liu, Sun, Li 2024 — scDesign3
*scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics.* Nature Biotechnology 2024;42:247–252 (published online 11 May 2023). doi:10.1038/s41587-023-01772-1. Peer reviewed (Brief Communication). R/Bioconductor package `scDesign3`. Source of arm M (PROTOCOLLO §5.2).

1. **Question:** Can one statistical model, fitted to real data, both generate realistic single-cell and spatial data (for benchmarking) and help interpret the real data (goodness-of-fit of clusters, trajectories, spatial locations)?
2. **Data:** Many public single-cell, multi-omics and spatial datasets; comparison with ZINB-WaVE, muscat and SPARSim ("three top-performing simulators for discrete cell types" in earlier benchmarks) and with a deep-learning simulator, scGAN, judged against left-out real cells.
3. **Main result:** For each gene, a GAMLSS model (a regression in which the mean, and possibly the dispersion and the zero inflation, depend on covariates) is fitted with a chosen distribution (Gaussian, Poisson, negative binomial, zero-inflated Poisson or zero-inflated NB). Cell type, pseudotime or spatial position enter as covariates, plus optional batch and condition. A copula (Gaussian, fast; or vine, slower) then adds the correlations between genes. Synthetic cells resembled left-out real cells better than the other simulators. One use shown in the paper: given single-cell and spot-level spatial data, scDesign3 generates spots with chosen cell-type proportions and was used to **benchmark deconvolution**, with results "consistent with a benchmark study" that CARD and RCTD outperformed SPOTlight in estimating proportions, while the three methods performed similarly in estimating each cell type's distribution across spots.
4. **Limit:** Benchmarks were run by the authors. Realism was judged on whole datasets; Crowell et al. 2023 showed that realism does not guarantee a faithful method ranking. Model choices (distribution, covariates, copula) are left to the user, and fitting gene by gene can be slow on large data.
5. **What we need from it:**
    - **What μ arm M takes:** the mean of each gene in each cell type from the fitted marginal model; the copula and the dispersion are not used, because our noise model is the same for all arms (§5.2).
    - **A point for the ADEMP review:** if the marginal model contains only cell type as a covariate, with a log link, the fitted mean of each cell type is equal to its average in the training data (the maximum-likelihood estimate in a regression with one categorical covariate is the group mean). In that case arm M would be almost identical to arm E. M differs from E only if the model includes other terms (for example library size, donor as batch, or a smoother), which must be decided and pre-registered. *(This reasoning is mine, not from the paper; to check in the scDesign3 vignette and with Luca.)*
    - **Precedent:** scDesign3 itself was used to simulate spots and rank deconvolution methods, which is the use that our Objective 2 tests.

---
*Benchmarks and candidate extra methods (2026-10-23).*

## Avila Cobos, Alquicira-Hernandez, Powell, Mestdagh, De Preter 2020 — bulk deconvolution benchmark
*Benchmarking of cell type deconvolution pipelines for transcriptomics data.* Nature Communications 2020;11:5650. doi:10.1038/s41467-020-19015-1. Peer reviewed.

1. **Question:** effect of transformation, normalization, marker selection, missing cell types and method on bulk deconvolution.
2. **Data:** pseudo-bulk mixtures from five scRNA-seq datasets; real PBMC bulk with known proportions.
3. **Main result:**
    - linear scale always best (log hurts);
    - best bulk methods: OLS, **nnls**, RLR, FARDEEP, CIBERSORT;
    - best single-cell-reference methods: **DWLS, MuSiC, SCDC**, as good as the best bulk ones;
    - a cell type missing from the reference worsens every method.
4. **Limit:** bulk, not spots; no spatial methods.
5. **What we need:** supports NNLS as a strong baseline and raw/linear counts as input; points to DWLS and MuSiC as candidates.

## Huuki-Myers, Montgomery, … Maynard, Collado-Torres 2025 — brain bulk deconvolution benchmark
*Benchmark of cellular deconvolution methods using a multi-assay dataset from postmortem human prefrontal cortex.* Genome Biology 2025;26:88. doi:10.1186/s13059-025-03552-3. Peer reviewed.

1. **Question:** which of six methods (DWLS, Bisque, MuSiC, BayesPrism, CIBERSORTx, hspe) best estimates cell-type proportions in human DLPFC?
2. **Data:** 22 DLPFC blocks; bulk RNA-seq (110 samples), snRNA-seq reference, real proportions from RNAScope/IF.
3. **Main result:**
    - **Bisque** (cor 0.54, RMSE 0.14) and **hspe** (0.53, 0.14) best, then CIBERSORTx;
    - MuSiC, BayesPrism and DWLS poor here (cor ≤ 0.05);
    - results depend on marker genes; Mean Ratio markers in R package `DeconvoBuddies`.
4. **Limit:** bulk, not spots; one tissue; rankings differ from Avila Cobos (DWLS and MuSiC good there, poor here).
5. **What we need:**
    - same region as our brain data (PsychAD DLPFC); `DeconvoBuddies` Mean Ratio could supply marker genes (e.g. SPOTlight `mgs`);
    - Bisque and hspe: not tested on spatial spots in any benchmark found; hspe is a v0.1 tarball (not on CRAN) → not proposed;
    - shows that benchmarks disagree on rankings, which is our thesis question.

## Spatial benchmarks checked online (for the method choice)
- **Li et al. 2022**, Nature Methods, doi:10.1038/s41592-022-01480-9: 32 simulated + real datasets; best deconvolution = **cell2location, SpatialDWLS, RCTD**.
- **Li et al. 2023**, Nature Communications 14:1548, doi:10.1038/s41467-023-37168-7 (a different paper from Li 2022): 18 methods, 50 real and simulated datasets; abstract: "CARD, Cell2location, and Tangram are the best methods for conducting the cellular deconvolution task" (checked 2026-10-26, as Luca asked). Used to state that CARD is excluded for design reasons, not for performance (DECISIONS 2026-10-26).
- **Sang-aram et al. 2024 (Spotless)**, eLife, doi:10.7554/eLife.88431: full entry below.
- **Conclusion (2026-10-23), now decided (DECISIONS 2026-10-26):** methods = RCTD, cell2location, spatialDWLS, MuSiC (second baseline), NNLS; CARD and SPOTlight are no longer in the study.

## Sang-aram, Browaeys, Seurinck, Saeys 2024 — Spotless
*Spotless, a reproducible pipeline for benchmarking cell type deconvolution in spatial transcriptomics.* eLife 2023;12:RP88431 (Version of Record 24 May 2024). doi:10.7554/eLife.88431. Peer reviewed (eLife reviewed preprint). PDF: `literature/Spotless_elife-88431-v1.pdf`. Independent benchmark chosen with Luca (DECISIONS 2026-10-26).

> **⚠ Still to be read in full by Pierangela.** Entry drafted by Claude on 2026-10-26 from the abstract, Results (pp. 1–7) and Appendix 1; to check and correct.

1. **Question:** Which spatial deconvolution methods are accurate, stable when the reference changes, and scalable?
2. **Data:** 11 methods. 63 "silver standards" made by their simulator synthspot from 7 scRNA-seq datasets (4 mouse brain, mouse kidney, mouse melanoma, human squamous cell carcinoma) × 9 abundance patterns, 10 replicates each, ~750 spots per replicate; 3 "gold standards" from seqFISH+ and STARMap (cells summed in 55 µm circles); 2 case studies (liver, melanoma). "Half of the cells from each scRNA-seq dataset were used to generate the synthetic spots and the other half used as the reference for deconvolution. This split was stratified by cell type."
3. **Main result:** "cell2location and RCTD are the top-performing methods, but surprisingly, a simple regression model outperforms almost half of the dedicated spatial deconvolution methods". Silver standards: RCTD and cell2location first, then SpatialDWLS, stereoscope and MuSiC; "All other methods ranked worse than NNLS in at least one metric". Gold standards: cell2location, MuSiC and RCTD the only ones above NNLS in all three metrics. Stability: "Except for MuSiC", the better methods were also more stable to a change of reference; "both NNLS and MuSiC were consistently in the bottom three". Metrics: RMSE, AUPR, JSD; rankings from the median over 10 replicates.
4. **Limit:** The split into reference and spots is **by cell** (stratified by cell type), not by donor, so reference and spots can share donors. Mostly mouse tissue. Performance varied a lot between replicates in the patterns with a dominant cell type. CARD not included.
5. **What we need from it:**
    - closest design to ours (single-cell data → synthetic spots → known truth), and the source of our baselines: NNLS beats about half of the methods; MuSiC strong on accuracy but unstable to the reference;
    - our split by donor is stricter than their split by cell: a point to make in Methods/Discussion;
    - same error measures as ours (RMSE, JSD); they also use AUPR (presence/absence);
    - replicates: 10 per setting, rankings on the median; our n_sim instead comes from the Monte Carlo SE (Morris);
    - their Appendix 1 validates the synthetic data with distributions of pairwise Spearman correlations (spots and genes, log CPM) and with mean–dispersion, library size and zero fractions: a possible model for our Check B.

## Dong, Yuan 2021 — SpatialDWLS (in the study, DECISIONS 2026-10-26)
*SpatialDWLS: accurate deconvolution of spatial transcriptomic data.* Genome Biology 2021;22:145. doi:10.1186/s13059-021-02362-7. R, inside `Giotto` (`runDWLSDeconv`).

1. **Question:** adapt DWLS (dampened weighted least squares, from bulk) to spots.
2. **Data:** 71 simulated spots from mouse-cortex seqFISH+ cells (523 cells, 6 cell types, ~51.5 µm squares); Visium mouse brain and developing human heart (no ground truth). Compared with MuSiC, RCTD, SPOTlight, stereoscope.
3. **Main result:** step 1 clusters spots and tests which cell types are enriched per cluster; step 2 DWLS only with those cell types. More accurate and faster than the others in their test (MuSiC and RCTD closest). Top 3 in Li 2022; just below RCTD/cell2location in Spotless.
4. **Limit:** tiny own validation (71 spots, 6 types, authors' benchmark); needs a spot clustering (`cluster_column`) and a signature matrix (both our choices); heavy `Giotto` dependency.
5. **What we need:**
    - no coordinates needed: `createGiottoObject()` makes dummy locations (checked in `GiottoClass` source, 2026-10-23); clustering is on expression;
    - `runDWLSDeconv(..., expression_values = "normalized", cluster_column = "leiden_clus", sign_matrix, n_cell = 50, cutoff = 2)` (Giotto source, `suite` branch); `n_cell` = cells per spot → set to our level;
    - would replace CARD (needs coordinates).
    - **Update 2026-10-26:** we use the standalone R package `spatialDWLS` (drieslab, v0.1.1, published 2026-10-16; no Giotto; imports Matrix, MatrixGenerics, quadprog, stats). Its DESCRIPTION says it follows Dong & Yuan 2021 in the same two rounds. Main calls in the vignette: `dwls_signature(ref_normalised, sign_gene, cell_type)` and `dwls(spot_expr_log, signature, clusters, n_cell)`. It needs a cluster per spot ("any clustering of the spots works (e.g. Leiden on the expression)"). It is a new re-implementation, not the Giotto version tested in Li 2022 and Spotless: a limit to declare (ADEMP §11).

## Wang, Park, Susztak, Zhang, Li 2019 — MuSiC (in the study as second baseline, DECISIONS 2026-10-26)
*Bulk tissue cell type deconvolution with multi-subject single-cell expression reference.* Nature Communications 2019;10:380. doi:10.1038/s41467-018-08023-x. R package `MuSiC` (GitHub xuranw/MuSiC, v1.0.0).

1. **Question:** bulk deconvolution that uses a single-cell reference from several subjects.
2. **Data:** pancreas and kidney bulk + multi-subject scRNA-seq.
3. **Main result:** weighted NNLS; genes weighted down if they vary a lot between subjects in the reference. Top in Avila Cobos; in Spotless above half of the spatial methods and top 3 on gold standards.
4. **Limit:** bulk method; poor in Huuki-Myers (brain); low stability to the reference in Spotless.
5. **What we need:**
    - no coordinates (spots = columns of `bulk.mtx`);
    - `music_prop(bulk.mtx, sc.sce, clusters, samples, markers = NULL, ...)` (source, 2026-10-23): needs per-cell **donor labels** (`samples`) → our training split has them; with few training donors the between-donor weights are rough (the number of training donors per tissue is fixed in ADEMP §2.5);
    - linear (raw) counts; genes with zero mean in the spots are dropped.

*Answer to the OSCA check (2026-10-23):* the built-in spot filters differ between methods: RCTD and CARD drop spots with fewer than 100 counts, cell2location drops none. None of the three published pipelines applies an adaptive (OSCA-style) spot filter. Keeping ADEMP §7.1 (no extra filter; spots a method drops count as failures) therefore matches what the published pipelines do. Whether to reopen the OSCA proposal is my decision.
