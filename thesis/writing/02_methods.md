# 2. Methods

<!--
Draft v2.0 (2026-11-01, for review). Same content as v1.0 (2026-10-31, from ADEMP v0.5 and DECISIONS up to 2026-10-30), reorganised and rewritten in the style of the Methods of a scientific paper, with thesis-level detail: standard section titles, factual reporting, short justifications with citations, equations and parameter tables. Internal references (DECISIONS rows, PROTOCOLLO sections, approval dates) are removed from the text; they go to the appendix on the protocol and its deviations (PROTOCOLLO §8 item 5). Hypotheses are in Section 2.8.1. Longer rationales from v1.0 are kept in the comment at the end of the file, for the Discussion.
Tense: present for the design, because the study has not been run yet; switch to the past tense when the results are in.
Placeholders: *[from the pilot]* = ADEMP [PILOT]; *[from the real data]* = ADEMP [REAL DATA]; *[pending Luca's approval]* = ADEMP [LUCA]/[OPEN]. [REF] marks a statement that still needs a checked source.
-->

## 2.1 Study design

The study is a simulation study planned and reported according to the ADEMP structure of Morris et al. (2019), which specifies aims, data-generating mechanisms, estimands, methods and performance measures. The protocol, including hypotheses, decision rules and the number of replicates, was fixed before any result was examined *[date of approval]*, and deviations from it are listed in Appendix B.

The study has two aims (Section 1.7). Aim 1 assesses fidelity, defined as the agreement between the cell-type mean expression profiles produced by each expression source and the profiles measured in donors that did not contribute to their construction. Aim 2 assesses the sensitivity of a deconvolution benchmark to the expression source. It asks whether the error and the ranking of five deconvolution methods change with the source, and which source yields the ranking closest to that obtained on pseudo-spots built from real cells of a held-out donor, referred to as the real anchor.

The principal factor is the source of the matrix μ of mean expression, with one row per gene and one column per cell type. Six sources, or arms, are compared (Section 2.4). All other components of the data-generating process are shared by the arms: the generator, the noise model, the spot compositions and the reference supplied to the methods. Three secondary factors are crossed with the principal factor in a fixed grid (Table 2.1). Five arms, three tissues and four combinations of cells per spot and depth give 60 data-generating mechanisms, or 72 if the generative arm enters the study (Section 2.4.4). Each of the twelve combinations of tissue, cells per spot and depth is matched by a real anchor. The hypotheses are stated in Section 2.8.1.

**Table 2.1.** Factors of the simulation study.

| Factor | Levels | Section |
|---|---|---|
| Expression source (μ) | E (empirical), P (parametric), M (model-based), G (generative), C− (negative control), E′ (independent comparator donor) | 2.4 |
| Tissue | Colorectal cancer (MMRp), thoracic lymph node, dorsolateral prefrontal cortex | 2.2 |
| Cells per spot | 2, 10 | 2.5.3 |
| Sequencing depth per cell | 175, 4,900 UMI | 2.5.3 |
| Deconvolution method | RCTD, cell2location, SpatialDWLS, MuSiC, NNLS | 2.6 |

## 2.2 Single-cell datasets

One single-cell RNA sequencing dataset is used per tissue. Within a tissue, the training, held-out and comparator donors are different donors of the same dataset, so that they share chemistry, processing pipeline and annotation, no integration across studies is required, and platform effects, which can otherwise dominate the biological signal in reference-based deconvolution (Cable et al. 2022), do not enter the comparisons. Datasets were searched in CZ CELLxGENE Discover, the Broad Single Cell Portal and the Gene Expression Omnibus. A dataset was eligible if it included at least ten donors and provided, for every cell, a cell-type label and a donor label assigned by the original authors, a requirement also adopted by Spotless for its single-cell references (Sang-aram et al. 2024). The cell-type label used throughout is the authors' annotation rather than the harmonised Cell Ontology term, and no cell is re-annotated. For colon, only healthy tissue or colorectal cancer was eligible. The three datasets are summarised in Table 2.2.

**Table 2.2.** Single-cell datasets.

| Tissue | Dataset | Technology | Subset used | Donors | Cells | Cell-type field |
|---|---|---|---|---|---|---|
| Brain | PsychAD (Fullard et al. 2025) | 10x Genomics 3′ v3, single-nucleus | Dorsolateral prefrontal cortex; disease status "normal"; age ≥ 18 years; donors of unknown age excluded | 208 | *[from `02_split.R`]* | `subclass` |
| Colon | Pelka et al. 2021 | 10x Genomics 3′ v2 | Tumours with proficient mismatch repair (MMRp) | 20 | 89,242 | `ClusterMidway` |
| Lymphoid tissue | Wells et al. 2025 | 10x Genomics 5′ v2 | Thoracic lymph nodes | 21 | 185,437 | `author_cell_type` |

Donor and cell counts after the subsetting, before donor quality control (Section 2.3.1), from the CELLxGENE cell metadata (accessed 27 October 2026) and the partition script (Section 2.9).

For brain, PsychAD provides single-nucleus profiles of healthy adult cortex. Single-nucleus data capture only part of the cellular mRNA but still resolve closely related cell types (Bakken et al. 2018). For colon, the dataset of Pelka et al. (2021) is restricted to MMRp tumours, because MMRp and mismatch-repair-deficient tumours differ both in composition and in the expression of the same cell types, and the original authors analysed them separately "to better capture differences between the two immunologically disparate tumor types" (Pelka et al. 2021). For lymphoid tissue, the atlas of Wells et al. (2025) is restricted to one anatomical site, the thoracic lymph node. This dataset contains immune cells only, so lymphoid pseudo-spots are mixtures of immune populations without stromal or endothelial cells. The three datasets differ in chemistry (3′ cells, 3′ nuclei and 5′ cells), and all comparisons are therefore made within tissue.

## 2.3 Data preprocessing

### 2.3.1 Donor quality control and partition

Donors with a low outlying number of cells are removed from each dataset before the partition. A donor is removed if its number of cells, on the logarithmic scale, is below the median minus three median absolute deviations, the threshold used by `scater::isOutlier` for the quality control of cells (Amezquita et al. 2020) applied here to donors. The rule removes seven brain donors, with 61 to 770 nuclei each and 3,335 nuclei in total, and no donor in colon or lymphoid tissue.

The data are partitioned by donor, never by cell. Among the donors that pass quality control, four are drawn at random (seed 20261029) to form the held-out group. The first donor drawn is the held-out donor of Aim 2, whose cells build the real anchor, and the second is the comparator donor of arm E′. All remaining donors form the training split (Table 2.3). In colon, all samples of a donor are assigned to the same split. The training split supplies arms E and M, the calibration of arm P, the library sizes and dispersions of the noise model and the reference given to every method. An automated test verifies that no donor identifier occurs in more than one split. If a step cannot process the full training split of brain, which contains about one million nuclei, cells, never donors, are subsampled within that split with a fixed seed.

**Table 2.3.** Donor partition and cell-type sets per tissue.

| | Brain | Colon | Lymphoid tissue |
|---|---|---|---|
| Donors removed by quality control | 7 | 0 | 0 |
| Training donors | 197 | 16 | 17 |
| Held-out group (in order of draw) | Donor_207, Donor_680, Donor_307, Donor_957 | C104, C105, C134, C162 | 583B, 591C, D543, D512 |
| Cell types with ≥ 50 training cells (sources and reference) | 27 | 18 | 30 |
| Cell types with ≥ 50 cells in the first held-out donor (Aim 1) | 20 | 8 | 9 |
| Cell types with ≥ 50 cells in training, held-out and comparator donor (spots) | 11 | 6 | 7 |

Source: partition script (Section 2.9). The first donor of the held-out group is the held-out donor of Aim 2; the second is the comparator donor of arm E′.

### 2.3.2 Cell-type sets

A cell type is used in a role only if it has at least 50 cells there, the threshold used by Crowell et al. (2023) to filter groups of their real references because insufficient numbers of cells interfere with the estimation of simulation parameters. Three sets of cell types follow in each tissue (Table 2.3). The expression sources and the reference of the methods use all types with at least 50 training cells. Aim 1 evaluates the types with at least 50 cells in the target held-out donor, and the remaining types are reported as not evaluated. The pseudo-spots of every arm and of the anchor contain only the types with at least 50 cells in training, in the held-out donor and in the comparator donor, so that all arms, E′ included, share the same compositions. The reference therefore contains more types than the spots, as a real atlas usually contains more types than a tissue section. Reference types absent from the spots have a true proportion of zero, and any proportion assigned to them counts as error. Types present in the spots but absent from the reference do not occur by construction; the omission of a cell type from the reference is known to increase deconvolution error substantially (Avila Cobos et al. 2020).

### 2.3.3 Gene and cell filters

Genes are retained if they have a count greater than one in at least ten training cells. The gene list is determined on the training split only and applied unchanged to all splits and arms, so that all matrices μ share the same genes. Cells with fewer than 100 detected genes are removed *[cell filter pending Luca's approval]*. Downloaded files are never modified, and all filtered data are written to separate files.

## 2.4 Expression sources

Each arm produces a matrix μ of dimension G × K, where G is the number of retained genes and K the number of cell types. Each column is a relative expression profile, that is the expected share of the transcriptome of a cell of type k taken by each gene. The six arms differ in the information from which μ is obtained (Table 2.4).

**Table 2.4.** Expression sources compared in the study.

| Arm | How μ is obtained | Data used | Random across replicates | Role in the study |
|---|---|---|---|---|
| E, empirical | Mean profile of the real training cells of each type | Training donors | No | Practice of using a real reference as source of profiles |
| P, parametric | Gene means drawn from a gamma distribution; marker genes multiplied by a fold-change | None, apart from the ranges matched to arm E | Yes | Fully synthetic profiles; maximal control, minimal realism |
| M, model-based | Cells simulated by scDesign3 fitted to the training cells (Aim 1: mean of the fitted model) | Training donors | Yes | Model-based simulator fitted to real data |
| G, generative | Mean profile predicted by GEM-1 from cell-type and tissue terms | Public archives used to train the model | Yes | Generative model trained on public data |
| C−, negative control | Profiles of arm E with gene labels permuted within each type | Training donors | Yes | Profiles without biological content; expected to fail |
| E′, comparator donor | Mean profile of the real cells of one other donor | One donor outside training | No | Reference for donor-to-donor variability |

The real anchor is not an arm: its pseudo-spots sum the real cells of the held-out donor (Section 2.5.4), and the rankings of all arms are compared with its ranking.

### 2.4.1 Empirical arm (E)

For cell type k, with C_k the set of its training cells and y_gc the count of gene g in cell c,

$$\mu^{E}_{gk} = \frac{\sum_{c \in C_k} y_{gc}}{\sum_{c \in C_k} L_c}, \qquad L_c = \sum_{g'} y_{g'c},$$

where the library size L_c is summed over all genes before filtering. The columns of μ^E restricted to the retained genes therefore sum to slightly less than one, and the mean count μ_gk L reproduces the real counts of the retained genes. This pooled estimate is the maximum-likelihood rate of a Poisson model with log L_c as offset and has the "profile times library size" structure used by RCTD (Cable et al. 2022) and muscat (Crowell et al. 2020). The same estimator is applied to the cells of the comparator donor (arm E′) and of each held-out donor (the targets of Aim 1).

### 2.4.2 Parametric arm (P)

Arm P generates μ de novo, without reference to real cells, following the approach of Splatter (Zappia et al. 2017). Baseline gene means are drawn from a gamma distribution, and for each cell type a set of marker genes is multiplied by a fold-change; the columns are then normalised to sum to one. The shape and rate of the gamma distribution, the number of markers per type and the fold-change are matched to the range of means and mean–variance relationship of arm E *[values from the real data; a priori values or calibration on arm E pending Luca's approval]*.

### 2.4.3 Model-based arm (M)

Arm M uses scDesign3 (Song et al. 2024) fitted to the training split, with cell type and donor as covariates. One model is fitted per tissue, since the training split excludes all donors of the held-out group; in brain, cells are subsampled with a fixed seed if the fit does not scale to the full split. Unlike the other arms, M is not reduced to a mean profile in Aim 2. A negative binomial model with cell type as covariate and library size as offset estimates the group mean almost exactly, so μ^M would coincide with μ^E and the arm would discard what distinguishes scDesign3 from the empirical estimator. The pseudo-spots of arm M are therefore built from cells generated by scDesign3, which carry its own marginal noise and gene–gene correlation, while compositions, cell types per spot and thinning are those of the other arms (Section 2.5). In each replicate the cells are generated for one training donor drawn at random, as the anchor comes from a single donor, and the donor covariate keeps between-donor variability out of the simulated noise, as for φ_g (Section 2.5.1). The difference between arms M and E thus measures the effect of the noise model at equal training data. In Aim 1, the profile of M is the mean of the fitted model over the training donors, computed once, and is expected to be close to μ^E by construction *[settings of marginal, copula and library size to be fixed in the pilot; donor covariate, one donor per replicate and Aim 1 profile pending Luca's approval]*.

### 2.4.4 Generative arm (G)

Arm G uses the mean profiles returned by the single-cell model of GEM-1 (gem-1-sc, Synthesize Bio) in its mean-estimation mode, queried with public metadata only, namely Cell Ontology terms for cell types and UBERON terms for tissues. The single-cell model was trained on CZ CELLxGENE Discover as of July 2024 (Koytiger et al. 2025). A dataset is therefore eligible for arm G if it is absent from the CELLxGENE Census of 1 July 2024, a condition met by all three datasets. The version of CELLxGENE underlying the model currently served by the provider's API is not known. Arm G enters the study only if written permission to publish, access to the API within the available quota and coverage of the tissues and cell types by the model's vocabulary are all obtained; this decision is recorded before any result of the arm is examined *[status of arm G]*. Every query is logged with its date, model identifier, full query and the MD5 checksum of the response.

### 2.4.5 Negative control (C−)

Arm C− takes μ^E and permutes the gene labels independently within each cell type. The distribution of expression values of each type is preserved and its biological content is removed. A new permutation is drawn in every replicate.

### 2.4.6 Independent comparator donor (E′)

Arm E′ applies the estimator of Section 2.4.1 to the cells of the comparator donor. No parameter is tuned on this donor. E′ provides a reference for donor-to-donor variability, against which the other arms are read, and is not an upper bound. Arms M and G estimate population-average profiles from many donors and may legitimately exceed E′. The attainable maximum of fidelity is estimated separately by a noise ceiling (Section 2.7.1).

## 2.5 Simulation of pseudo-spots

### 2.5.1 Noise model

The count of gene g in a simulated cell c of type k is drawn as

$$y_{gc} \sim \mathrm{NB}\left(m_{gc} = \mu_{gk} L_c,\; \phi_g\right), \qquad \mathrm{Var}(y_{gc}) = m_{gc} + \phi_g m_{gc}^2,$$

where μ_gk is taken from the arm, L_c is drawn at random from the library sizes of the training cells of type k, and φ_g is a gene-specific dispersion. The noise model is the same for all arms except M (Section 2.4.3), so that only μ differs between them. The negative binomial distribution describes droplet-based data with unique molecular identifiers without zero inflation (Svensson 2020), and a gene-specific dispersion is also assumed by muscat (Crowell et al. 2020) and by the model of cell2location (Kleshchevnikov et al. 2022). Empirical library sizes follow muscat, in which simulated cells take the library sizes of reference cells of their group, rather than a parametric distribution such as the log-normal of Splatter (Zappia et al. 2017).

The dispersions φ_g are estimated on the training split with glmGamPoi (Ahlmann-Eltze & Huber 2020), with design `~ cell_type + donor` and log L_c as offset. The donor term keeps between-donor variability out of φ_g. Without it, the simulated arms would be more variable than the anchor, which is built from a single donor, for a reason unrelated to μ. Per-gene maximum-likelihood estimates are used without shrinkage (`overdispersion_shrinkage = FALSE`), with the default Cox–Reid adjustment. In every tissue the fit uses 20,000 training cells sampled with a fixed seed and stratified by cell type in proportion to each type's share, because the brain training split does not fit in memory. Genes are fitted in parallel blocks, which gives estimates identical to a single fit. The adequacy of the sample is tested in colon, where all training cells can be fitted. The sample is accepted if the Pearson correlation of log φ_g between the sample and all cells is at least 0.99 across genes and the median ratio of the two estimates lies between 0.95 and 1.05; otherwise the sample size is raised to 50,000 cells in all tissues.

### 2.5.2 Spot composition

Each pseudo-spot dataset contains S spots *[S from the pilot; 500 in the pilot, the default of synthspot]*. The composition of spot s over the K_s cell types of the spots is drawn as

$$\boldsymbol{\pi}_s \sim \mathrm{Dirichlet}(\boldsymbol{\alpha}), \qquad \alpha_k = K_s \, p_k,$$

where p_k is the share of type k among the training cells of the spot types. With equal shares this reduces to the symmetric α = 1 used by the spatial-decomposition simulator of Open Problems [REF]. The n cells of the spot are then assigned to types by a multinomial draw with probabilities π_s, and the true composition is the fraction of cells of each type, π_sk = n_sk / n. Proportions taken from the training split make the spots resemble the tissue and limit the reuse of rare held-out cells in the anchor; they are taken from training rather than from the held-out donor so that the test donor does not shape its own test. Rare types consequently appear in few spots.

### 2.5.3 Cells per spot and sequencing depth

The number of cells per spot is fixed at n = 2 or n = 10, the bounds of the synthspot simulator used by Spotless (Sang-aram et al. 2024). Two cells represent an 8 µm bin of Visium HD, where bins typically mix two cell types (Oliveira et al. 2025), and ten cells a standard Visium spot, which contains 10 to 20 cells (Elosua-Bayes et al. 2021).

Sequencing depth is defined per cell, so that it remains crossed with the number of cells. The target depth of a spot is n times the depth per cell, with levels of 175 and 4,900 UMI per cell. The low level is derived from the mean UMI count per 8 µm bin under tissue in two public Visium HD human colorectal cancer samples of 10x Genomics, 207 and 488 UMI (Space Ranger metrics of `Visium_HD_Human_Colon_Cancer_P1` and `Visium_HD_Human_Colon_Cancer`), averaged and divided by two cells. The high level is derived from the median of 49,282 UMI per spot in the 10x Genomics sample `CytAssist_11mm_FFPE_Human_Colorectal_Cancer`, divided by ten cells [REF: 10x Genomics dataset pages]. The division by the number of cells assumes equal UMI per cell. The resulting targets are given in Table 2.5.

A spot whose total count exceeds its target is downsampled to exactly the target, by drawing the retained UMIs at random without replacement from all UMIs of the spot. Every UMI has the same probability of being retained, whichever cell and gene it comes from, so the expected count of each gene in each cell is reduced by the same factor, and the relative contribution of the cells to the spot, and of the genes within each cell, is preserved in expectation. A cell rich in RNA therefore continues to dominate the signal of a spot in proportion to its transcript content, as it does in tissue. Thinning represents the loss of transcripts between a cell and the spot that contains it. A spot covers only the part of each cell that lies within the section and within its capture area, and capture on the slide recovers a smaller fraction of the transcripts than droplet-based single-cell protocols [REF], so that the sum of whole-cell counts exceeds the depth of a real spot. Uniform thinning assumes that this loss affects all cells of a spot equally, whereas in tissue a spot may contain most of one cell and a small part of another; because the pseudo-spots have no geometry, this source of variation is not represented, a limitation considered in the Discussion.

A spot below its target is left unchanged, following synthspot. At the high level the realised depth of a spot is therefore the smaller of its target and the summed counts of its cells, and it varies with composition, since spots made of cell types with low RNA content stay below the target, as they would in tissue. At the low level every spot reaches its target. Thinning is applied identically to the simulated arms and to the anchor. Exact thinning replaced binomial thinning before the main run, because binomial thinning matches the target only in expectation and made the pre-registered depth check fail by chance in the pilot. For each level, arm and the anchor, the number of spots that reach the target, the number that stay below it and the median depth before thinning are reported *[from the pilot]*, to show whether the depth of the anchor differs from that of the arms. The lowest target, 350 UMI, is above the minimum of 100 UMI below which RCTD discards a spot.

**Table 2.5.** Target depth per spot (UMI).

| Cells per spot | 175 UMI per cell | 4,900 UMI per cell |
|---|---|---|
| 2 | 350 | 9,800 |
| 10 | 1,750 | 49,000 |

### 2.5.4 Real anchor

The anchor consists of pseudo-spots built by summing the counts of real cells of the held-out donor of Aim 2. For each spot, the n_sk cells of type k are drawn with replacement from the held-out cells of that type, using the same compositions as the simulated arms, and their counts are summed and thinned to the target depth. The true composition is therefore known while the expression is entirely real. Each anchor spot is built within a single sample. Because the held-out donors contain a few thousand cells, for example 2,433 in lymphoid donor 583B and 2,886 in colon donor C104, held-out cells are reused across spots and replicates, as in other pseudo-spot simulators; the expected reuse per cell is reported *[from the pilot]*.

### 2.5.5 Common random numbers

Within a replicate, all arms and the anchor share the same compositions and the same assignment of cell types to the cells of each spot, and only μ, and therefore the counts, differs, together with the noise model in arm M. Common random numbers induce a positive correlation between the results of different arms within a replicate and reduce the Monte Carlo variance of between-arm comparisons (Morris et al. 2019). Every method is applied to the same simulated datasets.

**Table 2.6.** Parameters of the generator.

| Parameter | Value | Source |
|---|---|---|
| μ_gk | Arm-specific (Section 2.4) | — |
| Library size L_c | Drawn from the training cells of type k | Crowell et al. 2020 |
| Dispersion φ_g | glmGamPoi, `~ cell_type + donor`, 20,000 training cells | Ahlmann-Eltze & Huber 2020 |
| Dirichlet α_k | K_s p_k, p_k = training share of type k | Open Problems [REF] |
| Cells per spot n | 2, 10 | synthspot; Oliveira et al. 2025; Elosua-Bayes et al. 2021 |
| Depth per cell | 175, 4,900 UMI | 10x Genomics public datasets |
| Spots per dataset S | *[from the pilot]* | synthspot default 500 in the pilot |

## 2.6 Deconvolution methods

Five methods are compared (Table 2.7). They were selected for their formulations and for their performance in independent benchmarks rather than for their popularity, which need not reflect performance (Crowell et al. 2023). RCTD and cell2location were the two best methods overall in Spotless (Sang-aram et al. 2024), and SpatialDWLS, cell2location and RCTD the three best in the benchmark of Li et al. (2022). The three dedicated methods represent three formulations: a Poisson likelihood with a correction for platform effects (RCTD), a Bayesian negative binomial model of absolute cell abundance (cell2location) and dampened weighted least squares (SpatialDWLS). Two baselines complete the set: non-negative least squares (NNLS) and MuSiC, a method developed for bulk tissue. CARD (Ma & Zhou 2022) is not included because it exploits the spatial correlation of composition between neighbouring spots, which pseudo-spots without positions do not have.

All methods receive, in every arm, the same reference, consisting of at most 1,000 training cells per cell type drawn at random with a fixed seed; types with fewer cells keep all of them. Methods are run with default parameters or the authors' recommended settings, and departures from the defaults are listed in Table 2.7. Built-in filters are left at their defaults, since they are part of the pipeline a user runs. The true number of cells per spot is supplied to the methods that accept it, which constitutes oracle information unavailable in real data.

RCTD is run with the R package spacexr in full mode, recommended by Cable et al. (2022) when spots may contain three or more cell types; the package default is doublet mode. Placeholder coordinates are supplied, since positions are not used in full mode. RCTD discards spots with fewer than 100 UMI and requires integer counts and at least 25 reference cells per type.

cell2location (Kleshchevnikov et al. 2022) is run in Python on a GPU with raw counts. The expected number of cells per location is set to the true n of each condition. The hyperparameter `detection_alpha` is set to 200 in the main analysis, the value recommended by the authors when detection sensitivity varies little across locations and used by Spotless for synthetic spots from the same technology as the reference, and to 20 in a sensitivity analysis *[pending Luca's approval]*.

SpatialDWLS (Dong & Yuan 2021) is run with the standalone R package spatialDWLS, which implements the method without the Giotto framework. Spots are clustered by k-means on principal components, candidate cell types are selected for each cluster by enrichment of marker genes, and each spot is then solved by dampened weighted least squares restricted to those types. Signature and marker genes are derived from the training reference, and the number of cells per spot is set to the true n *[clustering settings from the pilot]*.

MuSiC (Wang et al. 2019) weights genes by the consistency of their expression across the subjects of a multi-subject reference and receives the donor label of every reference cell. NNLS regresses each spot on the mean count per cell of each reference type under a non-negativity constraint, and the estimated coefficients are normalised to proportions.

**Table 2.7.** Deconvolution methods.

| Method | Implementation (version) | Model | Settings that differ from the defaults or depend on the condition | Built-in filter |
|---|---|---|---|---|
| RCTD | spacexr (2.2.1) | Poisson likelihood with platform-effect correction | Full mode; placeholder coordinates | Spots < 100 UMI; types < 25 cells |
| cell2location | cell2location *[version]* | Bayesian negative binomial | `N_cells_per_location` = n; `detection_alpha` = 200 (20 in sensitivity) | Genes in the reference |
| SpatialDWLS | spatialDWLS (0.1.1) | Dampened weighted least squares | `n_cell` = n; clustering *[from the pilot]* | *[to be verified]* |
| MuSiC | MuSiC (1.0.0) | Weighted non-negative least squares | Donor labels of the reference | Genes with zero mean in the spots |
| NNLS | nnls (1.6) | Non-negative least squares | Coefficients normalised to sum to one | None |

## 2.7 Performance measures

### 2.7.1 Fidelity (Aim 1)

The estimand of Aim 1 is the profile μ_·k of each cell type, compared with the profile of the same type in a held-out donor, μ^H_·k, estimated as in Section 2.4.1. Fidelity is computed for each arm, tissue and cell type rather than over all types jointly, because types differ in size and complexity (Crowell et al. 2023). The primary measure is the Pearson correlation of the log-profiles over the retained genes,

$$r_k = \mathrm{cor}\left(\log(\mu_{\cdot k} + \varepsilon),\; \log(\mu^{H}_{\cdot k} + \varepsilon)\right),$$

with the Spearman correlation as a sensitivity analysis *[pseudocount ε pending Luca's approval]*. Two secondary measures are the overlap of the top-k marker genes of each type and the separability of the cell types *[k and separability metric pending Luca's approval]*. Fidelity is computed against each of the four donors of the held-out group, and arm E′ is computed for all six pairs of donors within the group, which gives a distribution of donor-to-donor fidelity.

Because μ^H is itself estimated from a finite number of cells, a correlation of one cannot be reached even by a perfect source. The attainable maximum is estimated by a split-half noise ceiling, as in representational similarity analysis (Nili et al. 2014) [REF: verify]. For each tissue and cell type, the cells of the held-out donor are divided at random into two halves of equal size, the extra cell of an odd count being assigned to one half at random, and the correlation r_half between the two half-profiles is computed with the same coefficient and genes. The reliability of the full profile is obtained with the Spearman–Brown formula [REF], and its square root is the noise ceiling:

$$\rho_{\max} = \sqrt{\frac{2\,r_{\mathrm{half}}}{1 + r_{\mathrm{half}}}}.$$

The split is repeated B times and the results are averaged *[B and seed from the pilot]*.

### 2.7.2 Deconvolution error (Aim 2)

The estimand of Aim 2 is the proportion π_sk of cell type k in spot s. For method m in replicate i, with S spots and K_r cell types in the reference,

$$\mathrm{RMSE}_i = \sqrt{\frac{1}{S K_r}\sum_{s=1}^{S}\sum_{k=1}^{K_r}\left(\hat\pi_{sk} - \pi_{sk}\right)^2}, \qquad \mathrm{JSD}_i = \frac{1}{S}\sum_{s=1}^{S} \mathrm{JSD}\left(\hat{\boldsymbol{\pi}}_s,\, \boldsymbol{\pi}_s\right),$$

where π_sk = 0 for reference types absent from the spots and JSD is the Jensen–Shannon divergence. Each measure is averaged over the n_sim replicates and reported with its Monte Carlo standard error, MCSE = SD_i / √n_sim (Morris et al. 2019). Errors are computed for every combination of method, arm and condition and for the anchor. RMSE is the primary measure and JSD the secondary measure.

### 2.7.3 Rankings and their agreement

Within each arm and condition, the five methods are ranked by mean RMSE and, as a secondary analysis, by mean JSD. Agreement between the ranking of an arm and that of the anchor is measured with Kendall's τ, and agreement among the rankings of all arms with Kendall's W. Their MCSEs are estimated by bootstrap over replicates *[number of bootstrap samples from the pilot]*. With five methods, τ takes eleven values and a single exchange of two adjacent methods lowers it from 1 to 0.8; the pairwise differences in mean RMSE between methods are therefore reported with their MCSEs, so that exchanges within Monte Carlo error can be identified.

## 2.8 Statistical analysis

### 2.8.1 Hypotheses

Three hypotheses were stated before the study. H1: the negative control C− has the lowest fidelity in every tissue, and its ranking of methods differs from that of the anchor. H2: arm E′ quantifies the fidelity and the agreement with the anchor ranking attained by a second real donor, and provides the scale on which the other arms are interpreted. H3, the principal hypothesis: at least one synthetic arm (P, M or G) yields a ranking of methods that differs from the anchor ranking.

### 2.8.2 Decision rules

A ranking is considered to differ from the anchor ranking, or from the ranking of another arm, if at least one pairwise difference in mean RMSE between two methods changes sign and the difference exceeds twice its MCSE. The arm closest to the anchor is the arm with the highest τ, ties being broken by the mean absolute difference in RMSE from the anchor *[pending Luca's approval]*.

### 2.8.3 Number of replicates

The number of replicates is chosen to reach a target precision of the key performance measure (Morris et al. 2019):

$$n_{\mathrm{sim}} = \left(\frac{\mathrm{SD}}{\mathrm{MCSE}_{\mathrm{target}}}\right)^2,$$

where SD is the standard deviation of the per-replicate measure *[key measure, mean RMSE per method or between-method difference, and target MCSE pending Luca's approval]*. The SD is estimated in a pilot with ten replicates, so its relative standard error is about 1/√(2 × 9), or 24%. The upper limit of its 95% confidence interval is therefore used, approximately 1.83 times the estimate with nine degrees of freedom, and the largest value across the cells of the design is taken. The total computational cost is checked against the available resources *[from the pilot]*. cell2location, which runs on a shared GPU, may receive a separate number of replicates, fixed before any result is examined.

Arms E, E′ and M have no Monte Carlo error in Aim 1 once the partition is fixed, since their μ is not random (the profile of M is computed once from the fitted model); MCSEs of fidelity are reported for arms P, C− and G *[pending Luca's approval]*.

### 2.8.4 Pilot study

A pilot study, whose results are not part of the findings, provides the SD for n_sim, the computing time per run, the types of failure, the number of spots, the reuse of anchor cells and the check of seed handling (Section 2.9). It is run in colon, lymphoid tissue and brain, with arms E and C−, the anchor, the four methods that run on a CPU (RCTD, SpatialDWLS, MuSiC and NNLS) and ten replicates for each arm and condition. The pilot does not use the held-out group. In each tissue, a pilot held-out donor and a pilot comparator donor are drawn from the training donors (seed 20261030), and the pilot reference is built from the remaining training donors, so that all settings fixed by the pilot are chosen within training.

### 2.8.5 Failures

No quality filter is applied to the pseudo-spots, which are built from filtered cells; each method applies only its built-in filters. A failure is defined as a spot that a method discards or cannot process, a run that does not converge, or an estimated composition whose proportions do not sum to one within a tolerance *[tolerance pending Luca's approval]*. Failures are not discarded silently (Morris et al. 2019). They are counted for every combination of method, arm and condition and reported as a result. In the main analysis, methods are compared on the spots returned by all methods within a replicate; a method that fails on every spot of a replicate is excluded from this common set rather than emptying it. In a sensitivity analysis, failed spots receive the worst possible score and RMSE, τ and W are recomputed.

### 2.8.6 Control analyses

Three control analyses are planned. First, arm C− is expected to give low fidelity and rankings that differ from the anchor; a ranking unchanged under C− would indicate that the benchmark is insensitive to expression. Second, a synthetic arm whose fidelity exceeds the noise ceiling beyond Monte Carlo error is investigated as a possible error or leakage before the analysis continues; exceeding E′ alone is not taken as evidence of leakage. Third, data provenance is checked in all cases: no donor occurs in more than one split, and arm G is evaluated only on datasets absent from its training data.

## 2.9 Software, code and data availability

All analyses are performed in R 4.6.1, except cell2location, which runs in a dedicated Python environment *[environment and versions]*. R package versions are recorded with renv, and the environment is rebuilt from the lock file on the computing server. The main packages are glmGamPoi 1.24.0, spacexr 2.2.1, spatialDWLS 0.1.1, MuSiC 1.0.0, nnls 1.6 and rhdf5 2.56.1 *[scDesign3 and remaining versions from `renv.lock`]*. Data paths are read from an environment variable, so the same scripts run unchanged on a personal computer and on the server.

The pipeline consists of numbered scripts in the `thesis/` folder of the repository `lucavd/2025.geo_spatialtrans`, branch `thesis/expression-sources`: `01_download.R` downloads the data and records them in a manifest with URL, size, MD5 checksum and date; `02_split.R` performs donor quality control and the partition; `03_reference.R` computes μ for arms E and E′, library sizes, dispersions, the method reference and the held-out cells; `06_pilot.R` runs the pilot study. Functions are in `R/` (generator in `simulate_spots.R`, method wrappers in `deconvolve.R`, performance measures in `pilot_metrics.R`), and automated tests in `tests/` cover the partition, the generator and the dispersion sample. Downloaded data are not under version control and are never modified; on the server they are downloaded again and their checksums compared with the manifest. Queries to GEM-1 and the checksums of the responses are logged in `synthesize_log/`.

A single seed is fixed for the main study *[seed]*, and the state of the random number generator is stored at the start of each replicate (Morris et al. 2019). Parallel processes use independent streams of the L'Ecuyer-CMRG generator derived from that seed. The state of the generator is compared before and after two runs of each method in the pilot, to detect packages that reset the seed internally. Simulated datasets, or their seeds and generator states, are stored so that any replicate can be regenerated, and every figure and table is produced by a script from the stored data.

## References added in this chapter

Already in the Introduction: Avila Cobos et al. 2020; Bakken et al. 2018; Cable et al. 2022; Crowell et al. 2023; Dong & Yuan 2021; Kleshchevnikov et al. 2022; Koytiger et al. 2025; Li et al. 2022; Ma & Zhou 2022; Morris et al. 2019; Oliveira et al. 2025; Sang-aram et al. 2024; Song et al. 2024; Svensson 2020; Wang et al. 2019; Zappia et al. 2017.

New in this chapter (to be added to a common reference list):

Ahlmann-Eltze C, Huber W. glmGamPoi: fitting Gamma-Poisson generalized linear models on single cell count data. *Bioinformatics* 2020. doi:10.1093/bioinformatics/btaa1009

Amezquita RA, et al. Orchestrating single-cell analysis with Bioconductor. *Nature Methods* 2020. doi:10.1038/s41592-019-0654-x

Crowell HL, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nature Communications* 2020. doi:10.1038/s41467-020-19894-4

Elosua-Bayes M, et al. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. *Nucleic Acids Research* 2021;49:e50. doi:10.1093/nar/gkab043

Fullard JF, et al. *Scientific Data* 2025. doi:10.1038/s41597-025-04687-5 (PsychAD data release; full reference to check).

Li et al. 2023. *Nature Communications*. doi:10.1038/s41467-023-37168-7 (full author list to add).

Moreno-Torres JG, et al. A unifying view on dataset shift in classification. *Pattern Recognition* 2012;45:521–530. doi:10.1016/j.patcog.2011.06.019

Nili H, Wingfield C, Walther A, Su L, Marslen-Wilson W, Kriegeskorte N. A toolbox for representational similarity analysis. *PLoS Computational Biology* 2014. [REF: verify DOI and pages before citing].

Pelka K, Hofree M, Chen JH, et al. Spatially organized multicellular immune hubs in human colorectal cancer. *Cell* 2021;184:4734–4752. doi:10.1016/j.cell.2021.08.003

Wells et al. 2025. *Nature Immunology*. doi:10.1038/s41590-025-02241-4 (full reference to add).

Spearman–Brown prophecy formula: [REF: original source (Spearman 1910; Brown 1910) or a textbook].

10x Genomics Visium chemistries (poly(A) capture vs CytAssist probes) and the public datasets `Visium_HD_Human_Colon_Cancer`, `Visium_HD_Human_Colon_Cancer_P1`, `CytAssist_11mm_FFPE_Human_Colorectal_Cancer`: [REF: product documentation and dataset pages, accessed 2026-10-23 and 2026-10-31].

Open Problems spatial decomposition task (`openproblems-bio/task_spatial_decomposition`) and synthspot (`saeyslab/synthspot`): [REF: cite as software with access date].

<!--
Material moved out of Methods in v2.0 (from v1.0), to reuse in the Discussion or in Appendix B. Not thesis text.
- Why one dataset per tissue: platform effects quote (Cable 2022 p. 517); chemistries of Visium (poly(A) 3′ vs CytAssist probes) do not match any reference exactly → absolute errors transfer only partly to real spatial data (Discussion, limits).
- Brain nuclei / lymphoid 5′ acceptable because all comparisons are within tissue; between-tissue patterns confound tissue with chemistry and training size (limits).
- Colon MMRp: MMRp vs MMRd differ in T-cell compartment and in 16 of 43 malignant programmes (Pelka 2021); pooling would average two states; MMRp chosen over MMRd for ≥10 donors and more types with ≥50 cells (CELLxGENE metadata 2026-10-27). Colon results describe MMRp tumours, not CRC in general.
- Lymphoid thoracic only: one homogeneous site, more literature for biological checks; lymph nodes at different sites drain different organs.
- Arm G and chemistry: if the query cannot specify 5′, G may be disadvantaged in lymphoid tissue for technical reasons.
- Partition by donor vs by cell (Spotless): a cell split ignores donor variability and is optimistic; E′ measures that variability.
- Random draw of the held-out donor: typical donor, not the best sampled one. No rotation of the held-out donor (Steyerberg & Harrell 2016) → results hold for the drawn donors (limit, reduced by the held-out group of 4).
- Closed-set rule: unseen test class = dataset shift (Moreno-Torres 2012); Avila Cobos 2020 median RMSE up ≥2-fold in several settings when a type is removed.
- E′ not a ceiling: M and G average many donors and can be closer to a new donor than one other donor; P has no such explanation, so P above E′ is examined further.
- Dirichlet α from training: with α = 1 rare types get as many cells as common ones and the anchor reuses rare held-out cells many times; consequence: methods tested less on rare types (Spotless rare-type AUPR not included, BACKLOG).
- Depth: uniform thinning = partial capture of each cell (section, capture area, lower capture than droplets), relative contributions kept in expectation; limit = equal loss for all cells of a spot, no geometry. High level = min(target, summed cell counts), varies with composition; colon pilot 37–38% of 2-cell spots and 8–10% of 10-cell spots below target at 4,900 UMI per cell, equal across arms and anchor (06_pilot.R, depth CSV). Recheck in lymphoid and brain.
- RCTD doublet vs full mode with 2-cell spots: full mode handles them and keeps one setting across conditions.
- cell2location detection_alpha 20 is the authors' recommendation for strong detection gradients in real Visium data, which pseudo-spots lack.
- SpatialDWLS standalone package differs from the Giotto version evaluated in earlier benchmarks.
- MuSiC in Spotless: accurate but among the least stable when the reference changed. NNLS-type regression outperformed almost half of the dedicated methods in Spotless.
- CARD exclusion is not a performance judgement: among the best in Li et al. 2023.
- Earlier benchmarks fixed replicates without a precision target (Spotless 10; CARD 5 per setting; Crowell 2023 one dataset per reference × simulator, repository HelenaLC/simulation-comparison) → Introduction/Discussion.
- No QC on pseudo-spots: adaptive filters target damaged cells and empty droplets, absent from spots built from filtered cells; low-count spots are the designed low-depth level.
- Pilot within training "as a validation split in machine learning": test donors are used once.
- Stage reports (pre-registration, check C, check B, counter-proof, adversarial review) and dated amendments → Appendix B / reproducibility appendix.
-->
