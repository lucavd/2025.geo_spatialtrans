# 2. Methods

<!--
Draft v2.3 (2026-11-11): aligned with ADEMP v1.4 and DECISIONS of 2026-11-11. Brain pilot values (n_sim 33, final n_sim 169, null-baseline formula, CPU cost per tissue); marker cap of P (10,000 cells per type failed the lymph node check, so 20,000; P values of all tissues to recompute); streamed VST checked against Seurat, HVGs of P from its simulated reference; rule for the number of draws of P and C− and for B per tissue. Still open: lymph node and brain P values after the cap, depth design (question to Luca), B and draws per tissue, cell2location timing, 10x dataset pages.
Draft v2.2 (2026-11-10): review pass against DECISIONS up to 2026-11-06 and ADEMP v1.3. Brain cell count; cell filters before the partition; "lymph node" throughout; per-type fold-changes of P; Open Problems source; pilot numbers with the c rule (colon, lymph node: depth below target, n_sim, null-baseline formula with c); RMSE relative to the null baseline in Section 2.7.2; noise ceiling cited for the concept only (Nili et al. 2014), Spearman–Brown, Seurat (Stuart et al. 2019), Leiden settings, Benjamini–Hochberg and Pullin references checked; lymph node P calibration in Table 2.6; package versions from renv.lock; colons inside sentences removed. Still open: brain values, final n_sim (169 with the c rule vs 217 with c = 1), rule for the high depth level, 10x dataset pages.
Draft v2.1 (2026-11-03): aligned with DECISIONS of 2026-11-02 and 2026-11-03 (arm M with scDesign3 cells as a secondary analysis, exact thinning, failed spots filled with the null baseline, n_sim on paired differences with δ = 0.02, cell-share estimand with the UMI-share sensitivity) and with Pierangela's rules of 2026-11-03 from the literature and the colon pilot (parameters of P after OSCA, scDesign3 defaults, S = 500, pseudocount, top-100 markers with Jaccard, silhouette width, noise-ceiling B by MCSE, bootstrap B = 1,000, closest-arm rule, tolerance 1e-5 for numerical residues, uv environment, study seed 42); 2026-11-04: source check (P markers, spatialDWLS Leiden clustering, reference cap) and silhouette and Jaccard settings; 2026-11-04 later: aligned with ADEMP v1.1 and PROTOCOLLO 1.4 (P as a de novo world with a simulated reference, oracle arm O, arm M removed, H1 and H3 rewritten, cell2location replicates from a GPU budget, RMSE relative to the null baseline, P markers with `pval.type = "some"`, φ of P genes by rank). Values that depend on the lymphoid and brain pilots stay marked.
Draft v2.0 (2026-11-01, for review). Same content as v1.0 (2026-10-31, from ADEMP v0.5 and DECISIONS up to 2026-10-30), reorganised and rewritten in the style of the Methods of a scientific paper, with thesis-level detail: standard section titles, factual reporting, short justifications with citations, equations and parameter tables. Internal references (DECISIONS rows, PROTOCOLLO sections, approval dates) are removed from the text; they go to the appendix on the protocol and its deviations (PROTOCOLLO §8 item 5). Hypotheses are in Section 2.8.1. Longer rationales from v1.0 are kept in the comment at the end of the file, for the Discussion.
Tense: present for the design, because the study has not been run yet; switch to the past tense when the results are in.
Placeholders: *[from the pilot]* or *[values per tissue]* = values still to be computed by a fixed rule (lymphoid and brain pilots, calibration scripts); no rule is left open. [REF] marks a statement that still needs a checked source.
-->

## 2.1 Study design

The study is a simulation study planned and reported according to the ADEMP structure of Morris et al. (2019), which specifies aims, data-generating mechanisms, estimands, methods and performance measures. The protocol, including hypotheses, decision rules and the number of replicates, was fixed before any result was examined *[date of approval]*, and deviations from it are listed in Appendix B.

The study has two aims (Section 1.7). Aim 1 assesses fidelity, defined as the agreement between the cell-type mean expression profiles produced by each expression source and the profiles measured in donors that did not contribute to their construction. Aim 2 assesses the sensitivity of a deconvolution benchmark to the expression source. It asks whether the error and the ranking of five deconvolution methods change with the source, and which source yields the ranking closest to that obtained on pseudo-spots built from real cells of a held-out donor, referred to as the real anchor.

The principal factor is the source of the matrix μ of mean expression, with one row per gene and one column per cell type. Six sources, or arms, are compared (Section 2.4). All other components of the data-generating process are shared by the arms, namely the generator, the noise model, the spot compositions and, except in arm P, the reference supplied to the methods (Section 2.4.2). Three secondary factors are crossed with the principal factor in a fixed grid (Table 2.1). Five arms (E, P, C−, E′ and O), three tissues and four combinations of cells per spot and depth give 60 data-generating mechanisms, or 72 if the generative arm enters the study (Section 2.4.4). Each of the twelve combinations of tissue, cells per spot and depth is matched by a real anchor. The hypotheses are stated in Section 2.8.1.

**Table 2.1.** Factors of the simulation study.

| Factor | Levels | Section |
|---|---|---|
| Expression source (μ) | E (empirical), P (parametric), G (generative), C− (negative control), E′ (independent comparator donor), O (oracle) | 2.4 |
| Tissue | Colorectal cancer (MMRp), thoracic lymph node, dorsolateral prefrontal cortex | 2.2 |
| Cells per spot | 2, 10 | 2.5.3 |
| Sequencing depth per cell | 175, 4,900 UMI | 2.5.3 |
| Deconvolution method | RCTD, cell2location, SpatialDWLS, MuSiC, NNLS | 2.6 |

## 2.2 Single-cell datasets

One single-cell RNA sequencing dataset is used per tissue. Within a tissue, the training, held-out and comparator donors are different donors of the same dataset, so that they share chemistry, processing pipeline and annotation, no integration across studies is required, and platform effects, which can otherwise dominate the biological signal in reference-based deconvolution (Cable et al. 2022), do not enter the comparisons. Datasets were searched in CZ CELLxGENE Discover, the Broad Single Cell Portal and the Gene Expression Omnibus. A dataset was eligible if it included at least ten donors and provided, for every cell, a cell-type label and a donor label assigned by the original authors, a requirement also adopted by Spotless for its single-cell references (Sang-aram et al. 2024). The cell-type label used throughout is the authors' annotation rather than the harmonised Cell Ontology term, and no cell is re-annotated. For colon, only healthy tissue or colorectal cancer was eligible. The three datasets are summarised in Table 2.2.

**Table 2.2.** Single-cell datasets.

| Tissue | Dataset | Technology | Subset used | Donors | Cells | Cell-type field |
|---|---|---|---|---|---|---|
| Brain | PsychAD (Fullard et al. 2025) | 10x Genomics 3′ v3, single-nucleus | Dorsolateral prefrontal cortex; disease status "normal"; age ≥ 18 years; donors of unknown age excluded | 208 | 920,515 | `subclass` |
| Colon | Pelka et al. 2021 | 10x Genomics 3′ v2 | Tumours with proficient mismatch repair (MMRp) | 20 | 89,242 | `ClusterMidway` |
| Lymph node | Wells et al. 2025 | 10x Genomics 5′ v2 | Thoracic lymph nodes | 21 | 185,437 | `author_cell_type` |

Donor and cell counts after the subsetting and the cell filters (Section 2.3.3), before donor quality control (Section 2.3.1), from the CELLxGENE cell metadata (accessed 27 October 2026) and the partition script (Section 2.9).

For brain, PsychAD provides single-nucleus profiles of healthy adult cortex. Single-nucleus data capture only part of the cellular mRNA but still resolve closely related cell types (Bakken et al. 2018). For colon, the dataset of Pelka et al. (2021) is restricted to MMRp tumours, because MMRp and mismatch-repair-deficient tumours differ both in composition and in the expression of the same cell types, and the original authors analysed them separately "to better capture differences between the two immunologically disparate tumor types" (Pelka et al. 2021). For the lymph node, the atlas of Wells et al. (2025) is restricted to one anatomical site, the thoracic lymph node. This dataset contains immune cells only, so lymph node pseudo-spots are mixtures of immune populations without stromal or endothelial cells. The three datasets differ in chemistry (3′ cells, 3′ nuclei and 5′ cells), and all comparisons are therefore made within tissue.

## 2.3 Data preprocessing

### 2.3.1 Donor quality control and partition

Donors with a low outlying number of cells are removed from each dataset before the partition. A donor is removed if its number of cells, on the logarithmic scale, is below the median minus three median absolute deviations, the threshold used by `scater::isOutlier` for the quality control of cells (Amezquita et al. 2020) applied here to donors. The rule removes seven brain donors, with 61 to 770 nuclei each and 3,335 nuclei in total, and no donor in colon or lymph node.

The data are partitioned by donor, never by cell. Among the donors that pass quality control, four are drawn at random (seed 20261029) to form the held-out group. The first donor drawn is the held-out donor of Aim 2, whose cells build the real anchor, and the second is the comparator donor of arm E′. All remaining donors form the training split (Table 2.3). In colon, all samples of a donor are assigned to the same split. The training split supplies arm E, the calibration of arm P, the library sizes and dispersions of the noise model and the reference given to every method. An automated test verifies that no donor identifier occurs in more than one split. If a step cannot process the full training split of brain, which contains about 900,000 nuclei, cells, never donors, are subsampled within that split with a fixed seed.

**Table 2.3.** Donor partition and cell-type sets per tissue.

| | Brain | Colon | Lymph node |
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

Before the partition, cells with fewer than 100 detected genes and cells without an author cell-type label are removed. After the partition, genes are retained if they have a count greater than one in at least ten training cells. The gene list is determined on the training split only and applied unchanged to all splits and arms, so that all matrices μ share the same genes. The gene and cell thresholds follow Crowell et al. (2023), who filtered their real references before estimating simulation parameters. Downloaded files are never modified, and all filtered data are written to separate files.

## 2.4 Expression sources

Each arm produces a matrix μ of dimension G × K, where G is the number of retained genes and K the number of cell types. Each column is a relative expression profile, that is the expected share of the transcriptome of a cell of type k taken by each gene. The arms differ in the information from which μ is obtained (Table 2.4).

**Table 2.4.** Expression sources compared in the study.

| Arm | How μ is obtained | Data used | Random across replicates | Role in the study |
|---|---|---|---|---|
| E, empirical | Mean profile of the real training cells of each type | Training donors | No | Practice of using a real reference as source of profiles |
| P, parametric | Gene means drawn from a gamma distribution; marker genes multiplied by a fold-change | Training donors, only to calibrate the distributions | Yes | Fully synthetic world, with a method reference simulated from μ^P; maximal control, minimal realism |
| G, generative | Mean profile predicted by GEM-1 from cell-type and tissue terms | Public archives used to train the model | Yes | Generative model trained on public data |
| C−, negative control | Profiles of arm E with gene labels permuted within each type | Training donors | Yes | Profiles without biological content; expected to fail |
| E′, comparator donor | Mean profile of the real cells of one other donor | One donor outside training | No | Reference for donor-to-donor variability |
| O, oracle | Mean profile of the real cells of the held-out donor itself | Held-out donor | No | Reference for the effect of the noise model (Aim 2 only) |

The real anchor is not an arm. Its pseudo-spots sum the real cells of the held-out donor (Section 2.5.4), and the rankings of all arms are compared with its ranking.

### 2.4.1 Empirical arm (E)

For cell type k, with C_k the set of its training cells and y_gc the count of gene g in cell c,

$$\mu^{E}_{gk} = \frac{\sum_{c \in C_k} y_{gc}}{\sum_{c \in C_k} L_c}, \qquad L_c = \sum_{g'} y_{g'c},$$

where the library size L_c is summed over all genes before filtering. The columns of μ^E restricted to the retained genes therefore sum to slightly less than one, and the mean count μ_gk L reproduces the real counts of the retained genes. This pooled estimate is the maximum-likelihood rate of a Poisson model with log L_c as offset and has the "profile times library size" structure used by RCTD (Cable et al. 2022) and muscat (Crowell et al. 2020). The same estimator is applied to the cells of the comparator donor (arm E′) and of each held-out donor (the targets of Aim 1).

### 2.4.2 Parametric arm (P)

Arm P generates μ de novo, following the approach of Splatter (Zappia et al. 2017). Baseline gene means are drawn from a gamma distribution, and for each cell type a set of marker genes is multiplied by a fold-change; the columns are then normalised to sum to one. No real profile enters μ^P; the training split is used only to calibrate the three distributions, so that the profiles of P cover the same range of expression as those of arm E. Marker genes of each type are identified on the training cells following Amezquita et al. (2020) with `scran::findMarkers`, which applies Welch t-tests between each pair of types, blocked on donor, testing up-regulation against a log2 fold-change threshold of one, an approach similar to TREAT (McCarthy & Smyth 2009). A gene is a marker of a type if it is up-regulated against at least half of the other types (`pval.type = "some"`) with a Benjamini–Hochberg false discovery rate (Benjamini & Hochberg 1995) below 0.05. The stricter requirement of up-regulation against every other type misses genes shared by related types (Amezquita et al. 2020); in the colon training data it left seven of the 18 types without markers, four of them types present in the pilot pseudo-spots, and these types would have had identical profiles in P. Because the cell-type labels were derived from the same data, these p-values serve to rank genes rather than as strict inference (Amezquita et al. 2020); the number of markers obtained at a false discovery rate of 0.01 is reported alongside. Because the brain training split (897,347 nuclei) does not fit in memory, marker detection uses at most 20,000 cells per type in every tissue, drawn at random within each type with a fixed seed; types with fewer cells enter whole, so that every type is represented. The gene filter, the baseline profile and the fold-changes are computed on all training cells. A first cap of 10,000 cells, the default number of cells per type retained in the reference of RCTD (Cable et al. 2022), was checked in the lymph node, which has four types above 10,000 cells and fits in memory, by the Jaccard index between the marker sets obtained with and without the cap, required in advance to be at least 0.9 for every type. The check failed in six of 30 types, with a Jaccard of 0.60 for central memory CD4 T cells, which had four markers in both runs, and between 0.75 and 0.86 in five other types, while the number of markers per type changed by at most eight. As specified in advance, the cap was therefore doubled to 20,000 in every tissue; at 20,000 the Jaccard index was at least 0.9 in 28 of 30 types, 0.75 for central memory CD4 T cells (three markers with the cap, four without) and 0.86 for MAIT cells. MAIT cells were below the cap and entered whole, so the remaining differences arise from the cells of the other types against which each type is tested. The baseline distribution is estimated with the procedure of Splatter on the pooled training profile of the genes that are not markers of any type. The gamma distribution is fitted by minimising the Cramér–von Mises distance to the means winsorized at their 10th and 90th percentiles, and genes whose log mean exceeds the median by more than two median absolute deviations are treated as expression outliers, whose share and log-normal factors relative to the median mean are estimated separately (Zappia et al. 2017). Splatter averages normalised cells, whereas we use the pooled profile on which the generator is defined; our implementation returns the same parameters as the Splatter functions on the same input. In each replicate, baseline means are drawn as in Splatter, with outlier genes set to the median mean multiplied by a log-normal factor. The fold-change of a training marker is the mean profile of its type divided by the average profile of the other types. Every type then receives as many markers as were found for it in the training data, drawn at random among the genes, and each marker is multiplied by a fold-change drawn with replacement from the fold-changes of the training markers of the same type, so that each type keeps the strength of its own markers. A gene drawn as a marker of two types receives both fold-changes, and the sharing of markers between related types is therefore not reproduced beyond chance. Because only the number and strength of the markers of each type are taken from the training data, arm P reproduces cell types defined by a few strong markers but not the many small shifts in expression that separate closely related states. In the lymph node, four of the seven cell types of the pseudo-spots are CD4 T-cell states with three to six markers each, and their separability in P may therefore be higher or lower than in the real cells; it is measured by the per-type silhouette width of Aim 1 (Section 2.7.1). A new μ^P is drawn in every replicate. The calibrated values are reported in Table 2.6 *[values per tissue from the calibration script; colon done]*.

Because the genes of P do not correspond to real genes, pseudo-spots generated from μ^P would resemble none of the cell types of a real reference, and arm P would reduce to a second negative control. Arm P is therefore built as a coherent synthetic world. In every replicate the reference supplied to the methods is simulated from the same μ^P, with one simulated cell for each cell of the real reference (Section 2.6), carrying its cell type, donor label and library size, and counts drawn from the common noise model (Section 2.5.1). Arm P thus asks whether a benchmark built entirely de novo ranks the methods as a benchmark on real cells does. The dispersion of a gene of P, which has no real counterpart, is that of the real gene with the same rank of mean expression, where genes of P are ranked by their mean baseline value across types and real genes by their pooled training mean. Overdispersion depends on the mean, as in the biological coefficient of variation of Splatter (Zappia et al. 2017), and matching by rank preserves both this relation and the real scatter of the dispersions among genes of similar mean.

### 2.4.3 Oracle arm (O)

Arm O applies the estimator of Section 2.4.1 to the cells of the held-out donor itself, the donor from which the real anchor is built, and generates pseudo-spots with the common noise model. It is a reference arm rather than a competing source. Its μ is the true mean profile of the anchor donor, so the difference between the anchor and O measures the effect of the simulated noise model, and the difference between O and the other arms measures the effect of the source of μ. Arm O is used in Aim 2 only, since in Aim 1 its profile coincides with the held-out profile by construction.

### 2.4.4 Generative arm (G)

Arm G uses the mean profiles returned by the single-cell model of GEM-1 (gem-1-sc, Synthesize Bio) in its mean-estimation mode, queried with public metadata only, namely Cell Ontology terms for cell types and UBERON terms for tissues. The single-cell model was trained on CZ CELLxGENE Discover as of July 2024 (Koytiger et al. 2025). A dataset is therefore eligible for arm G if it is absent from the CELLxGENE Census of 1 July 2024, a condition met by all three datasets. The version of CELLxGENE underlying the model currently served by the provider's API is not known. Arm G enters the study only if written permission to publish, access to the API within the available quota and coverage of the tissues and cell types by the model's vocabulary are all obtained; this decision is recorded before any result of the arm is examined *[status of arm G]*. Every query is logged with its date, model identifier, full query and the MD5 checksum of the response.

### 2.4.5 Negative control (C−)

Arm C− takes μ^E and permutes the gene labels independently within each cell type. The distribution of expression values of each type is preserved and its biological content is removed. A new permutation is drawn in every replicate.

### 2.4.6 Independent comparator donor (E′)

Arm E′ applies the estimator of Section 2.4.1 to the cells of the comparator donor. No parameter is tuned on this donor. E′ provides a reference for donor-to-donor variability, against which the other arms are read, and is not an upper bound. Arm G estimates population-average profiles from many donors and may legitimately exceed E′. The attainable maximum of fidelity is estimated separately by a noise ceiling (Section 2.7.1).

## 2.5 Simulation of pseudo-spots

### 2.5.1 Noise model

The count of gene g in a simulated cell c of type k is drawn as

$$y_{gc} \sim \mathrm{NB}\left(m_{gc} = \mu_{gk} L_c,\; \phi_g\right), \qquad \mathrm{Var}(y_{gc}) = m_{gc} + \phi_g m_{gc}^2,$$

where μ_gk is taken from the arm, L_c is drawn at random from the library sizes of the training cells of type k, and φ_g is a gene-specific dispersion. The noise model is the same for all arms, so that only μ differs between them. The negative binomial distribution describes droplet-based data with unique molecular identifiers without zero inflation (Svensson 2020), and a gene-specific dispersion is also assumed by muscat (Crowell et al. 2020) and by the model of cell2location (Kleshchevnikov et al. 2022). Empirical library sizes follow muscat, in which simulated cells take the library sizes of reference cells of their group, rather than a parametric distribution such as the log-normal of Splatter (Zappia et al. 2017).

The dispersions φ_g are estimated on the training split with glmGamPoi (Ahlmann-Eltze & Huber 2020), with design `~ cell_type + donor` and log L_c as offset. The donor term keeps between-donor variability out of φ_g. Without it, the simulated arms would be more variable than the anchor, which is built from a single donor, for a reason unrelated to μ. Per-gene maximum-likelihood estimates are used without shrinkage (`overdispersion_shrinkage = FALSE`), with the default Cox–Reid adjustment. In every tissue the fit uses 20,000 training cells sampled with a fixed seed and stratified by cell type in proportion to each type's share, because the brain training split does not fit in memory. Genes are fitted in parallel blocks, which gives estimates identical to a single fit. The adequacy of the sample is tested in colon, where all training cells can be fitted. The sample is accepted if the Pearson correlation of log φ_g between the sample and all cells is at least 0.99 across genes and the median ratio of the two estimates lies between 0.95 and 1.05; otherwise the sample size is raised to 50,000 cells in all tissues.

### 2.5.2 Spot composition

Each pseudo-spot dataset contains S = 500 spots, the default of the synthspot simulator (Sang-aram et al. 2024), in every tissue and condition. The composition of spot s over the K_s cell types of the spots is drawn as

$$\boldsymbol{\pi}_s \sim \mathrm{Dirichlet}(\boldsymbol{\alpha}), \qquad \alpha_k = c \, p_k,$$

where p_k is the share of type k among the training cells of the spot types and c is the concentration, which sets how mixed the spots are. A small c gives spots dominated by one type, a large c gives spots close to the tissue average. With c = K_s and equal shares the model would reduce to the symmetric α = 1 used by the spatial-decomposition simulator of Open Problems (repository `openproblems-bio/task_spatial_decomposition`, which accepts a vector α). We instead set c from real tissue, separately for each tissue and number of cells per spot. In the public 10x Genomics Visium HD sections Human Colon Cancer, Human Lymph Node and Mouse Brain [REF: 10x dataset pages], nuclei with a cell type assigned by marker enrichment were grouped into sets of 2 and 10 neighbours within regions of 1 mm². In each region c was taken as the value at which the variance of the compositions of these sets under a Dirichlet-multinomial model equals the observed variance. The value used is the median over the regions of a tissue, five per tissue and ten in the colon, where tumour and stromal regions of the same section are pooled because a tumour sample contains both. The median is preferred to the mean because in the brain at 10 nuclei a single region gives c = 110. The resulting values at 2 and 10 cells are 0.68 and 1.05 in colon, 1.79 and 2.05 in lymph node, and 3.37 and 21.8 in brain. In every region c is lower at 2 nuclei than at 10, because two neighbouring nuclei are more often of the same type than ten. These sections were typed with a coarser catalogue (three to six types per tissue) than the author labels used here, in the lymph node about half of the nuclei could not be typed, and the brain sections are from mouse, whereas our brain data are human. The n cells of the spot are then assigned to types by a multinomial draw with probabilities π_s, and the true composition is the fraction of cells of each type, π_sk = n_sk / n. Proportions taken from the training split make the spots resemble the tissue and limit the reuse of rare held-out cells in the anchor; they are taken from training rather than from the held-out donor so that the test donor does not shape its own test. Because the α_k are proportional to p_k, whatever the value of c, the expected proportion of type k in a spot equals its share in the tissue, E[π_sk] = p_k. Rare types consequently appear in few spots.

### 2.5.3 Cells per spot and sequencing depth

The number of cells per spot is fixed at n = 2 or n = 10, the bounds of the synthspot simulator used by Spotless (Sang-aram et al. 2024). Two cells represent an 8 µm bin of Visium HD, where bins typically mix two cell types (Oliveira et al. 2025), and ten cells a standard Visium spot, which contains 10 to 20 cells (Elosua-Bayes et al. 2021).

Sequencing depth is defined per cell, so that it remains crossed with the number of cells. The target depth of a spot is n times the depth per cell, with levels of 175 and 4,900 UMI per cell. The low level is derived from the mean UMI count per 8 µm bin under tissue in two public Visium HD human colorectal cancer samples of 10x Genomics, 207 and 488 UMI (Space Ranger metrics of `Visium_HD_Human_Colon_Cancer_P1` and `Visium_HD_Human_Colon_Cancer`), averaged and divided by two cells. The high level is derived from the median of 49,282 UMI per spot in the 10x Genomics sample `CytAssist_11mm_FFPE_Human_Colorectal_Cancer`, divided by ten cells [REF: 10x Genomics dataset pages]. The division by the number of cells assumes equal UMI per cell. The resulting targets are given in Table 2.5.

A spot whose total count exceeds its target is downsampled to exactly the target, by drawing the retained UMIs at random without replacement from all UMIs of the spot. Every UMI has the same probability of being retained, whichever cell and gene it comes from, so the expected count of each gene in each cell is reduced by the same factor, and the relative contribution of the cells to the spot, and of the genes within each cell, is preserved in expectation. A cell rich in RNA therefore continues to dominate the signal of a spot in proportion to its transcript content, as it does in tissue. Thinning represents the loss of transcripts between a cell and the spot that contains it. A spot covers only the part of each cell that lies within the section and within its capture area, and capture on the slide recovers a smaller fraction of the transcripts than droplet-based single-cell protocols [REF], so that the sum of whole-cell counts exceeds the depth of a real spot. Uniform thinning assumes that this loss affects all cells of a spot equally, whereas in tissue a spot may contain most of one cell and a small part of another; because the pseudo-spots have no geometry, this source of variation is not represented, a limitation considered in the Discussion.

A spot below its target is left unchanged, following synthspot. At the high level the realised depth of a spot is therefore the smaller of its target and the summed counts of its cells, and it varies with composition, since spots made of cell types with low RNA content stay below the target, as they would in tissue. At the low level every spot reaches its target. Thinning is applied identically to the simulated arms and to the anchor. Exact thinning replaced binomial thinning before the main run, because binomial thinning matches the target only in expectation and made the pre-registered depth check fail by chance in the pilot. In the pilots, at 4,900 UMI per cell, 40–42% of the two-cell spots and 18–22% of the ten-cell spots stayed below the target in colon, in similar proportions in the arms and in the anchor. In the lymph node the shares were higher and differed between the anchor (44% of the two-cell and 33% of the ten-cell spots) and the simulated arms (58–59% and 53–54%), because many lymph node cells have small libraries *[brain from the pilot; rule for the high level under review]*. For each level, arm and the anchor, the number of spots that reach the target, the number that stay below it and the median depth before thinning are reported, to show whether the depth of the anchor differs from that of the arms. The lowest target, 350 UMI, is above the minimum of 100 UMI below which RCTD discards a spot.

**Table 2.5.** Target depth per spot (UMI).

| Cells per spot | 175 UMI per cell | 4,900 UMI per cell |
|---|---|---|
| 2 | 350 | 9,800 |
| 10 | 1,750 | 49,000 |

### 2.5.4 Real anchor

The anchor consists of pseudo-spots built by summing the counts of real cells of the held-out donor of Aim 2. For each spot, the n_sk cells of type k are drawn with replacement from the held-out cells of that type, using the same compositions as the simulated arms, and their counts are summed and thinned to the target depth. The true composition is therefore known while the expression is entirely real. Each anchor spot is built within a single sample. Because the held-out donors contain a few thousand cells, for example 2,433 in lymph node donor 583B and 2,886 in colon donor C104, held-out cells are reused across spots and replicates, as in other pseudo-spot simulators; the expected reuse per cell is reported *[from the pilot]*.

### 2.5.5 Common random numbers

Within a replicate, all arms and the anchor share the same compositions and the same assignment of cell types to the cells of each spot, and only μ, and therefore the counts, differs. Common random numbers induce a positive correlation between the results of different arms within a replicate and reduce the Monte Carlo variance of between-arm comparisons (Morris et al. 2019). Every method is applied to the same simulated datasets.

**Table 2.6.** Parameters of the generator.

| Parameter | Value | Source |
|---|---|---|
| μ_gk | Arm-specific (Section 2.4) | — |
| Library size L_c | Drawn from the training cells of type k | Crowell et al. 2020 |
| Dispersion φ_g | glmGamPoi, `~ cell_type + donor`, 20,000 training cells | Ahlmann-Eltze & Huber 2020 |
| Dirichlet α_k | c p_k, p_k = training share of type k | Open Problems `task_spatial_decomposition`; synthspot |
| Concentration c (2 / 10 cells) | colon 0.68 / 1.05, lymph node 1.79 / 2.05, brain 3.37 / 21.8; median over regions of Visium HD neighbourhoods of real nuclei | 10x Visium HD public datasets [REF: dataset pages] |
| Cells per spot n | 2, 10 | synthspot; Oliveira et al. 2025; Elosua-Bayes et al. 2021 |
| Depth per cell | 175, 4,900 UMI | 10x Genomics public datasets |
| Spots per dataset S | 500 | synthspot default (Sang-aram et al. 2024) |
| Arm P: gamma shape and rate | Splatter procedure (Cramér–von Mises on winsorized means) on the pooled training profile of non-marker genes; counts per 10,000; shape 0.588 and rate 4.07 in colon, 0.546 and 2.48 in lymph node, *[brain]* | Zappia et al. 2017 |
| Arm P: expression outliers | Share and log-normal factor of genes with log mean > median + 2 MAD; share 1.7%, location 4.65, scale 1.05 in colon, 1.25%, 5.05 and 0.90 in lymph node, *[brain]* | Zappia et al. 2017 |
| Arm P: dispersion of a gene | φ of the real gene with the same rank of mean expression | Zappia et al. 2017 (mean-dependent dispersion) |
| Arm P: markers per type and fold-changes | As found in training by `findMarkers` (`lfc = 1`, `pval.type = "some"`, FDR < 0.05) on at most 20,000 cells per type; 638 marker genes in colon (16–144 per type) and 341 in lymph node (3–74), fold-change median 5.5 and 5.0; *[brain]* | Amezquita et al. 2020 |

## 2.6 Deconvolution methods

Five methods are compared (Table 2.7). They were selected for their formulations and for their performance in independent benchmarks rather than for their popularity, which need not reflect performance (Crowell et al. 2023). RCTD and cell2location were the two best methods overall in Spotless (Sang-aram et al. 2024), and SpatialDWLS, cell2location and RCTD the three best in the benchmark of Li et al. (2022). The three dedicated methods represent three formulations, namely a Poisson likelihood with a correction for platform effects (RCTD), a Bayesian negative binomial model of absolute cell abundance (cell2location) and dampened weighted least squares (SpatialDWLS). Two baselines complete the set, non-negative least squares (NNLS) and MuSiC, a method developed for bulk tissue. CARD (Ma & Zhou 2022) is not included because it exploits the spatial correlation of composition between neighbouring spots, which pseudo-spots without positions do not have.

All methods receive, in every arm except P (Section 2.4.2), the same reference, consisting of at most 1,000 training cells per cell type drawn at random with a fixed seed; types with fewer cells keep all of them. RCTD and Spotless allow up to 10,000 cells per type; the smaller cap limits the cost of the methods that work on dense matrices and loses little information, since in the colon pilot the profile of each capped type correlated with the profile of all its training cells with a Pearson coefficient of at least 0.998. The same check is applied in every tissue, and the cap is raised to 10,000 if any type falls below 0.99. Methods are run with default parameters or the authors' recommended settings, and departures from the defaults are listed in Table 2.7. Built-in filters are left at their defaults, since they are part of the pipeline a user runs. No parameter is tuned on the simulated data. The true number of cells per spot is supplied to the methods that accept it, which constitutes oracle information unavailable in real data.

With the normalisation recommended by their authors, the methods do not all estimate the same quantity. RCTD defines its weights as "the proportion of the contribution of cell type k to pixel i" (Cable et al. 2022), that is the share of the transcripts of a spot contributed by each type, and SpatialDWLS builds its signature on expression normalised per cell and rescales its solution to one, so both estimate shares of UMI. MuSiC and NNLS, whose profiles are mean counts per cell, estimate shares of cells. The two quantities differ when cell types differ in RNA content; in the training cells of the colon pilot, the median library size was 12,453 UMI in epithelial tumour cells and 2,146–2,298 UMI in CD8 and CD4 T cells. All methods are scored against the share of cells (Section 2.7.2), and the gap is reported as a property of a method whose output is read as a cell proportion. The gap depends only on the library sizes of the cells, which are the same in every arm, so it does not bias the comparison between sources of μ.

RCTD is run with the R package spacexr in full mode, recommended by Cable et al. (2022) when spots may contain three or more cell types; the package default is doublet mode. Placeholder coordinates are supplied, since positions are not used in full mode. RCTD discards spots with fewer than 100 UMI and requires integer counts and at least 25 reference cells per type.

cell2location (Kleshchevnikov et al. 2022) is run in Python on a GPU with raw counts. The expected number of cells per location is set to the true n of each condition. The hyperparameter `detection_alpha` is set to 200 in the main analysis, the value recommended by the authors when detection sensitivity varies little across locations and used by Spotless for synthetic spots from the same technology as the reference, and to 20 in a sensitivity analysis.

SpatialDWLS (Dong & Yuan 2021) is run with the standalone R package spatialDWLS, which implements the method without the Giotto framework. Spots are clustered, candidate cell types are selected for each cluster by enrichment of marker genes, and each spot is then solved by dampened weighted least squares restricted to those types. The signature is built from the training reference on normalised, non-logged expression, with the 100 genes of highest log2 fold-change per type, as in Spotless (Sang-aram et al. 2024). Spots are log-normalised and clustered with the settings of the Spotless pipeline (repository of Sang-aram et al. 2024, `subworkflows/deconvolution/spatialdwls/script_nf.R`), in which a shared-nearest-neighbour graph with four neighbours on the first ten principal components of the marker genes is partitioned by the Leiden algorithm with resolution 0.4 and 1,000 iterations (Traag et al. 2019). The graph is built with bluster and the clustering with igraph rather than Giotto, and principal components are computed on the marker genes rather than on highly variable genes. The parameter `n_cell` is set to the true n, the expected number of cells per spot as defined by the authors. The method keeps a cell type for a cluster only if its proportion in some spot reaches 1/`n_cell`; with two cells per spot this threshold equals the true proportion of a type, so whole clusters may return no estimate. This behaviour is the same in the Giotto implementation and in the standalone package, and the resulting failures are reported as a property of the method at Visium HD resolution, where 8 µm bins hold one or two cells. A sensitivity analysis uses `n_cell` = 50, the default of Giotto.

MuSiC (Wang et al. 2019) weights genes by the consistency of their expression across the subjects of a multi-subject reference and receives the donor label of every reference cell. NNLS regresses each spot on the mean count per cell of each reference type under a non-negativity constraint, and the estimated coefficients are normalised to proportions.

**Table 2.7.** Deconvolution methods.

| Method | Implementation (version) | Model | Settings that differ from the defaults or depend on the condition | Built-in filter |
|---|---|---|---|---|
| RCTD | spacexr (2.2.1) | Poisson likelihood with platform-effect correction | Full mode; placeholder coordinates | Spots < 100 UMI; types < 25 cells |
| cell2location | cell2location *[version]* | Bayesian negative binomial | `N_cells_per_location` = n; `detection_alpha` = 200 (20 in sensitivity) | Genes in the reference |
| SpatialDWLS | spatialDWLS (0.1.1) | Dampened weighted least squares | `n_cell` = n (50 in sensitivity); top 100 markers per type; Leiden clustering (k = 4 neighbours, 10 principal components, resolution 0.4) | None on spots; cell types below 1/`n_cell` in a cluster |
| MuSiC | MuSiC (1.0.0) | Weighted non-negative least squares | Donor labels of the reference | Genes with zero mean in the spots |
| NNLS | nnls (1.6) | Non-negative least squares | Coefficients normalised to sum to one | None |

## 2.7 Performance measures

### 2.7.1 Fidelity (Aim 1)

The estimand of Aim 1 is the profile μ_·k of each cell type, compared with the profile of the same type in a held-out donor, μ^H_·k, estimated as in Section 2.4.1. Fidelity is computed for each arm, tissue and cell type rather than over all types jointly, because types differ in size and complexity (Crowell et al. 2023). The primary measure is the Pearson correlation of the log-profiles over the retained genes,

$$r_k = \mathrm{cor}\left(\log\left(1 + 10^4\,\mu_{\cdot k}\right),\; \log\left(1 + 10^4\,\mu^{H}_{\cdot k}\right)\right),$$

that is the correlation of the log-transformed profiles expressed as counts per 10,000 with a pseudocount of one, the default log-normalisation of Seurat (`NormalizeData`, scale factor 10,000; Stuart et al. 2019). The Spearman correlation is used as a sensitivity analysis.

Two secondary measures follow. The first is the overlap between the marker genes of a type in the source and in the held-out donor. For each type, the 100 genes with the highest log2 fold-change against the mean of the other types, computed on the same log-normalised scale as the correlation, are selected in each of the two sets of profiles, and their overlap is measured by the Jaccard index, the size of the intersection divided by the size of the union. One hundred markers per type is the number used for marker-based deconvolution in Spotless (Sang-aram et al. 2024); scDesign3 used 50 (Song et al. 2024), and in a benchmark of marker selection methods the number of genes, between 5 and 40, had little effect on the ranking of the methods (Pullin & McCarthy 2024). Because no single value is established, the Jaccard index is also reported as a function of the number of markers, for 10, 20, 50, 100 and 200 genes, as in the evaluation of GEM-1 (Koytiger et al. 2025). The second is the separability of the cell types, measured as in Crowell et al. (2023) by the silhouette width (Rousseeuw 1987) with cell type as the group label, computed on Euclidean distances in the space of the principal components. Cells are generated from each source with the common noise model (Section 2.5.1); each type receives as many cells as it has in the held-out donor. Their mean silhouette width per type is compared with that of the real cells of the held-out donor. Counts are log-normalised with scuttle, a shifted logarithm that performs as well as or better than more sophisticated transformations (Ahlmann-Eltze & Huber 2023). The 2,000 most highly variable genes are selected once on the raw counts of the training split with the variance-stabilising method of Seurat (`FindVariableFeatures`, `selection.method = "vst"`, with its default of 2,000 genes; Stuart et al. 2019) and used for every arm and for the held-out donor, so that the held-out data do not influence the choice; the first 50 principal components are computed with scater, as in Crowell et al. (2023). In a recent benchmark of feature selection methods, Seurat-VST obtained the highest overall ranking and 2,000 features consistently produced high scores (Zappia et al. 2025); this number lies within the range of 1,000 to 5,000 recommended by Luecken & Theis (2019). Crowell et al. (2023) selected 500 genes with scran, a method that did not rank among the best in that benchmark. Because the brain training split does not fit in memory, the variance-stabilising selection is computed from per-gene sums accumulated over blocks of cells; on the colon training cells this implementation selected the same 2,000 genes as Seurat 6.0.0, with standardised variances differing by less than 10⁻¹². The genes of arm P have no counterpart among these genes, so for P the 2,000 genes are selected with the same method on its simulated reference (Section 2.4.2), in each draw.

Fidelity is computed against each of the four donors of the held-out group, and arm E′ is computed for all six pairs of donors within the group, which gives a distribution of donor-to-donor fidelity.

Because μ^H is itself estimated from a finite number of cells, a correlation of one cannot be reached even by a perfect source. The attainable maximum is estimated by a noise ceiling, the expected performance of a perfect model given the noise in the data, a concept used in representational similarity analysis (Nili et al. 2014), here computed from split-half reliability. For each tissue and cell type, the cells of the held-out donor are divided at random into two halves of equal size, the extra cell of an odd count being assigned to one half at random, and the correlation r_half between the two half-profiles is computed with the same coefficient and genes. The reliability of the full profile is obtained with the Spearman–Brown formula (Spearman 1910; Brown 1910), and its square root is the noise ceiling:

$$\rho_{\max} = \sqrt{\frac{2\,r_{\mathrm{half}}}{1 + r_{\mathrm{half}}}}.$$

The split is repeated B times and the ceilings are averaged. B is chosen so that the Monte Carlo standard error of the averaged ceiling, SD(ρ_max)/√B, is below 0.005 for every cell type, below the second decimal place at which fidelity is reported. As for the number of replicates (Section 2.8.3), the SD is taken as the upper limit of its 95% confidence interval from a pilot of ten splits, and B = ⌈(SD_upper/0.005)²⌉ with a minimum of ten, computed in each tissue *[B per tissue]*.

The genes of arm P do not correspond to real genes, so its gene-by-gene correlation and marker overlap with the held-out profiles are close to zero by construction. Arm P is therefore evaluated only on distributional properties, each compared with the held-out donor and with arm E, namely the separability of its cell types, measured by the silhouette width as above; the distribution of the gene means, summarised by the 10th, 50th, 90th, 99th and 99.9th percentiles of log(1 + 10⁴ μ) over genes and by the two-sample Kolmogorov–Smirnov distance; and the distribution of the fold-changes of the top 100 markers of each type, summarised in the same way.

### 2.7.2 Deconvolution error (Aim 2)

The estimand of Aim 2 is the proportion π_sk of cell type k in spot s, defined as a share of cells and used as the truth for every method. For method m in replicate i, with S spots and K_r cell types in the reference,

$$\mathrm{RMSE}_i = \sqrt{\frac{1}{S K_r}\sum_{s=1}^{S}\sum_{k=1}^{K_r}\left(\hat\pi_{sk} - \pi_{sk}\right)^2}, \qquad \mathrm{JSD}_i = \frac{1}{S}\sum_{s=1}^{S} \mathrm{JSD}\left(\hat{\boldsymbol{\pi}}_s,\, \boldsymbol{\pi}_s\right),$$

where π_sk = 0 for reference types absent from the spots and JSD is the Jensen–Shannon divergence. Each measure is averaged over the n_sim replicates and reported with its Monte Carlo standard error, MCSE = SD_i / √n_sim (Morris et al. 2019). Errors are computed for every combination of method, arm and condition and for the anchor. RMSE is the primary measure and JSD the secondary measure. Each mean RMSE is also reported relative to that of the null baseline, which assigns to every spot the share of each cell type among all training cells, the best guess without information, in the same arm and condition. A ratio below one indicates that the method adds information beyond the expected composition. The ratio is accompanied by the paired difference between the RMSE of the method and that of the null baseline within each replicate, with its MCSE, so that a ratio close to one is read against its Monte Carlo error. The failure rate of each method is a further secondary measure (Section 2.8.5).

For the sensitivity analysis on the scale of the method outputs (Section 2.6), the share of UMI of type k in spot s is computed from the library sizes of the cells of the spot,

$$\pi^{\mathrm{UMI}}_{sk} = \frac{\sum_{c \in s,\, k(c) = k} L_c}{\sum_{c \in s} L_c},$$

and RCTD and SpatialDWLS are scored against this truth, MuSiC and NNLS against the share of cells.

### 2.7.3 Rankings and their agreement

Within each arm and condition, the five methods are ranked by mean RMSE and, as a secondary analysis, by mean JSD. Agreement between the ranking of an arm and that of the anchor is measured with Kendall's τ, and agreement among the rankings of all arms with Kendall's W. Their MCSEs are estimated by a bootstrap over replicates with B = 1,000 resamples, more than the 200 usually sufficient for a standard error and enough for percentile intervals (Efron & Tibshirani 1993, chapters 6 and 19). With five methods, τ takes eleven values and a single exchange of two adjacent methods lowers it from 1 to 0.8; the pairwise differences in mean RMSE between methods are therefore reported with their MCSEs, so that exchanges within Monte Carlo error can be identified.

## 2.8 Statistical analysis

### 2.8.1 Hypotheses

Three hypotheses were stated before the study. Under H1, among the arms whose genes are real genes, the negative control C− has the lowest fidelity in every tissue, and its ranking of methods differs from that of the anchor. Arm P is excluded from this comparison, because its gene-by-gene agreement with the held-out profiles is null by construction; in Aim 1 it is evaluated on the distributions of gene means and fold-changes and on cell-type separability. Under H2, arm E′ quantifies the fidelity and the agreement with the anchor ranking attained by a second real donor, and provides the scale on which the other arms are interpreted. Under H3, the principal hypothesis, the ranking of methods obtained with arm P or arm G differs from the anchor ranking more than the rankings of the reference arms O and E′ do.

### 2.8.2 Decision rules

A pair of methods changes sign between an arm and the anchor if the difference in mean RMSE between the two methods has opposite signs in the arm and in the anchor and both differences exceed twice their MCSE. A ranking is considered to differ from the anchor ranking if at least one pair changes sign. For H3, sign changes against the anchor are counted per arm over method pairs, conditions and tissues; H3 is supported for P or G if its count exceeds the counts of both O and E′. With ten pairs, four conditions and three tissues, some sign changes are expected by chance, and the counts of the reference arms O and E′ provide this background. Among the source arms E, P and G, the arm closest to the anchor is the arm with the highest τ, ties being broken by the mean absolute difference in RMSE from the anchor.

### 2.8.3 Number of replicates

The number of replicates is chosen to reach a target precision of the key performance measure (Morris et al. 2019). Because rankings are decided by differences between methods, the key measure is the paired difference in RMSE between two methods within the same replicate, D_i = RMSE_i(m_1) − RMSE_i(m_2), which under common random numbers varies less than the RMSE of a single method. The smallest relevant difference between two methods is set to δ = 0.02, and the MCSE of D is required to satisfy 2 × MCSE ≤ δ/2, so that

$$n_{\mathrm{sim}} = \left(\frac{4\,\mathrm{SD}_{\mathrm{upper}}(D)}{\delta}\right)^2.$$

The SD is estimated in a pilot with ten replicates, so its relative standard error is about 1/√(2 × 9), or 24%. The upper limit of its 95% confidence interval, SD_upper, is therefore used, approximately 1.83 times the estimate with nine degrees of freedom, and the largest value across method pairs and cells of the design is taken. In the pilots, SD_upper(D) was 0.065 in colon and 0.019 in the lymph node, which give n_sim = 169 and 14; in colon all of the largest values involved SpatialDWLS. With δ = 0.02, 73% of the pairs of methods on arm E and on the anchor differed by at least δ in colon and 58% in the lymph node. In the brain it was 0.029 (n_sim = 33, with 69% of the pairs differing by at least δ). The rule is applied to the largest SD across the three tissues, which gives n_sim = 169, with a minimum of 72 replicates if the cost of cell2location allows. The total computational cost, n_sim × data-generating mechanisms × methods × time per run, is checked against the available resources. For the methods that run on a CPU, one replicate took 1.2, 2.5 and 11.0 CPU-hours in colon, lymph node and brain, with MuSiC accounting for 68% of the time; memory limits the brain to five parallel processes *[design of the depth levels and resulting run time]*. For cell2location the time is measured on the server *[cell2location timing]*. cell2location runs on a GPU with a budget of 240 hours. The time of one complete replicate, comprising all pseudo-datasets of the replicate and the training of the reference of arm P in each tissue, is measured first, and cell2location receives n_c2l = 240 h divided by this time replicates, at most the n_sim of the other methods and at least 20; below 20, it is run only at ten cells per spot. Paired differences involving cell2location are computed on its first n_c2l replicates *[n_c2l from the server timing]*. This number is fixed before any result is examined.

Arms E and E′ have no Monte Carlo error in Aim 1 once the partition is fixed, since their μ is not random; MCSEs of fidelity are reported for arms P and C−, and for arm G if repeated queries return different profiles. The number of draws of P and C− follows the rule of the noise ceiling (Section 2.7.1). The SD of each fidelity measure, per cell type and held-out donor, is estimated from ten draws, the upper limit of its 95% confidence interval is taken, and the number of draws is ⌈(SD_upper/0.005)²⌉, the largest value over types and measures, with a minimum of ten; it is computed in each tissue *[number of draws per tissue]*.

### 2.8.4 Pilot study

A pilot study, whose results are not part of the findings, provides the SD for n_sim, the computing time per run, the types of failure, the number of spots, the reuse of anchor cells and the check of seed handling (Section 2.9). It is run in colon, lymph node and brain, with arms E and C−, the anchor, the four methods that run on a CPU (RCTD, SpatialDWLS, MuSiC and NNLS) and ten replicates for each arm and condition. The pilot does not use the held-out group. In each tissue, a pilot held-out donor and a pilot comparator donor are drawn from the training donors (seed 20261030), and the pilot reference is built from the remaining training donors, so that all settings fixed by the pilot are chosen within training.

### 2.8.5 Failures

No quality filter is applied to the pseudo-spots, which are built from filtered cells; each method applies only its built-in filters. A failure is defined as a spot that a method discards or cannot process, a run that does not converge, or an estimated composition with missing values, a proportion below −10⁻⁵ or proportions whose sum differs from one by more than 10⁻⁵. Values between −10⁻⁵ and zero are treated as numerical zeros, set to zero, and the proportions are rescaled to one. The tolerance lies above the residues observed in the pilot, where RCTD returned two negative values out of about 1.08 million, none below −7 × 10⁻⁷, and above the rounding error expected from single-precision arithmetic, used by cell2location, over 30 cell types (about 3.6 × 10⁻⁶). Failures are not discarded silently (Morris et al. 2019). They are counted for every combination of method, arm and condition and reported as a result.

No method is excluded from a condition because of failures, and no threshold on the failure rate is applied. The failure rate depends on the arm, so it is an outcome rather than a filter, and a threshold would yield rankings with different sets of methods in different arms, which Kendall's τ could not compare. Excluding a method would also conceal poor performance linked to failure (Pawel et al. 2026). In the main analysis, every failed spot receives as its estimate the share of each reference cell type among all training cells, the best prediction available without information on the spot, used as a null baseline in Spotless (Sang-aram et al. 2024). Failures therefore count in the RMSE, and spots and methods are the same in every arm. The fill affects a method's error only to the extent that the method, where it does return an estimate, performs better than the null baseline. For this reason, the error of the null baseline applied to every spot is reported as a reference value, the RMSE of every method is also reported relative to it (the ratio of the two mean RMSEs in the same arm and condition, below one when the method adds information), and the failure rate is always read alongside the RMSE. Under the composition model of Section 2.5.2 the expected composition of every spot equals the training shares, so the null baseline is close to the best constant prediction, and its error is set by the spread of the true compositions around their mean. With concentration c and n cells per spot, the expected squared error for type k is p_k(1 − p_k)[1/(c + 1) + c/((c + 1)n)], which decreases with n and increases as c decreases. In the pilots this expression gave an RMSE of 0.182 at two cells and 0.149 at ten in colon, 0.137 and 0.105 in the lymph node, and 0.134 and 0.066 in the brain, in agreement with the observed error of the null baseline in every arm and condition. The failure rate, the fraction of spots on which a method fails within a replicate, is a secondary performance measure, averaged over replicates with its Monte Carlo standard error and compared with the rate on the anchor, which gives the level of failure expected on real cells.

Three sensitivity analyses are planned, in each of which RMSE, τ and W are recomputed. In the first, methods are compared on the spots returned by all methods within a replicate; a method that fails on every spot of a replicate is left out of this common set rather than emptying it. In the second, SpatialDWLS is run with `n_cell` = 50 (Section 2.6). In the third, each method is scored against the truth that its output represents, the share of UMI for RCTD and SpatialDWLS and the share of cells for MuSiC and NNLS (Section 2.7.2). Because n_sim is set on the main analysis, the third analysis is less precise *[SD_upper(D) against the UMI-share truth from the pilots]*.

### 2.8.6 Control analyses

Three control analyses are planned. First, arm C− is expected to give low fidelity and rankings that differ from the anchor; a ranking unchanged under C− would indicate that the benchmark is insensitive to expression. Second, a synthetic arm whose fidelity exceeds the noise ceiling beyond Monte Carlo error is investigated as a possible error or leakage before the analysis continues; exceeding E′ alone is not taken as evidence of leakage. Third, data provenance is checked in all cases, so that no donor occurs in more than one split and arm G is evaluated only on datasets absent from its training data.

## 2.9 Software, code and data availability

All analyses are performed in R 4.6.1 (4.6.2 on the computing server), except cell2location, which runs in a dedicated Python environment created with uv, whose package versions are fixed in a lock file kept under version control *[cell2location version]*. R package versions are recorded with renv, and the environment is rebuilt from the lock file on the computing server. The main packages are glmGamPoi 1.24.0, scran 1.40.0, scuttle 1.22.0, bluster 1.22.0, igraph 2.3.4, spacexr 2.2.1, spatialDWLS 0.1.1, MuSiC 1.0.0, nnls 1.6 and rhdf5 2.56.1 *[scater and Seurat versions once installed]*. Data paths are read from an environment variable, so the same scripts run unchanged on a personal computer and on the server.

The pipeline consists of numbered scripts in the `thesis/` folder of the repository `lucavd/2025.geo_spatialtrans`, branch `thesis/expression-sources`. `01_download.R` downloads the data and records them in a manifest with URL, size, MD5 checksum and date, `02_split.R` performs donor quality control and the partition, `03_reference.R` computes μ for arms E and E′, library sizes, dispersions, the method reference and the held-out cells, `04p_calibrate_p.R` calibrates arm P, `04_simulate.R` generates the pseudo-spots of all arms and the anchor, `05_deconvolve.R` runs the methods, and `06_pilot.R` and `06b_pilot_nsim.R` run the pilot study and compute n_sim. Functions are in `R/` (generator in `simulate_spots.R`, arm P in `arm_p.R`, method wrappers in `deconvolve.R`, performance measures in `pilot_metrics.R`), and automated tests in `tests/` cover the partition, the generator and the dispersion sample. Downloaded data are not under version control and are never modified; on the server they are downloaded again and their checksums compared with the manifest. Queries to GEM-1 and the checksums of the responses are logged in `synthesize_log/`.

All random steps use fixed seeds written in the code. The donor partition (seed 20261029), the pilot (20261030 and 20261101) and the sampling of the method reference and of the cells for φ_g (20261031) use the seeds fixed when those steps were run. All later steps, namely the main simulation, the calibration of arm P, the noise ceiling and the bootstrap, use the study seed 42, and the state of the random number generator is stored at the start of each replicate (Morris et al. 2019). Parallel processes use independent streams of the L'Ecuyer-CMRG generator derived from the study seed, with a fixed number of workers. The state of the generator was compared after the runs of each method in the colon pilot, and no method reset it *[lymph node and brain from the pilot]*. Simulated datasets, or their seeds and generator states, are stored so that any replicate can be regenerated, and every figure and table is produced by a script from the stored data.

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

Nili H, Wingfield C, Walther A, Su L, Marslen-Wilson W, Kriegeskorte N. A toolbox for representational similarity analysis. *PLoS Computational Biology* 2014;10:e1003553. doi:10.1371/journal.pcbi.1003553.

Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* 1995;57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x.

McCarthy DJ, Smyth GK. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 2009;25:765–771. doi:10.1093/bioinformatics/btp053

Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 2019;9:5233. doi:10.1038/s41598-019-41695-z.

Efron B, Tibshirani RJ. *An Introduction to the Bootstrap*. Chapman & Hall, 1993.

Pullin JM, McCarthy DJ. A comparison of marker gene selection methods for single-cell RNA sequencing data. *Genome Biology* 2024;25:56. doi:10.1186/s13059-024-03183-0.

Ahlmann-Eltze C, Huber W. Comparison of transformations for single-cell RNA-seq data. *Nature Methods* 2023. doi:10.1038/s41592-023-01814-1

Zappia L, Richter S, Ramírez-Suástegui C, et al. Feature selection methods affect the performance of scRNA-seq data integration and querying. *Nature Methods* 2025;22:834. doi:10.1038/s41592-025-02624-3

Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology* 2019;15:e8746. doi:10.15252/msb.20188746

Rousseeuw PJ. Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. *Journal of Computational and Applied Mathematics* 1987;20:53–65.

Pawel S, Bartoš F, Siepe BS, Lohmann A. Handling missingness, failures, and non-convergence in simulation studies: a review of current practices and recommendations. *The American Statistician* 2026;80:31–48. doi:10.1080/00031305.2025.2540002

Pelka K, Hofree M, Chen JH, et al. Spatially organized multicellular immune hubs in human colorectal cancer. *Cell* 2021;184:4734–4752. doi:10.1016/j.cell.2021.08.003

Wells et al. 2025. *Nature Immunology*. doi:10.1038/s41590-025-02241-4 (full reference to add).

Spearman C. Correlation calculated from faulty data. *British Journal of Psychology* 1910;3:271–295.

Brown W. Some experimental results in the correlation of mental abilities. *British Journal of Psychology* 1910;3:296–322.

Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM 3rd, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive integration of single-cell data. *Cell* 2019;177:1888–1902. doi:10.1016/j.cell.2019.05.031.

10x Genomics Visium chemistries (poly(A) capture vs CytAssist probes) and the public datasets `Visium_HD_Human_Colon_Cancer`, `Visium_HD_Human_Colon_Cancer_P1`, `CytAssist_11mm_FFPE_Human_Colorectal_Cancer`: [REF: product documentation and dataset pages, accessed 2026-10-23 and 2026-10-31].

Open Problems spatial decomposition task (`openproblems-bio/task_spatial_decomposition`) and synthspot (`saeyslab/synthspot`): [REF: cite as software with access date].

<!--
Material moved out of Methods in v2.0 (from v1.0), to reuse in the Discussion or in Appendix B. Not thesis text.
- Why one dataset per tissue: platform effects quote (Cable 2022 p. 517); chemistries of Visium (poly(A) 3′ vs CytAssist probes) do not match any reference exactly → absolute errors transfer only partly to real spatial data (Discussion, limits).
- Brain nuclei / lymphoid 5′ acceptable because all comparisons are within tissue; between-tissue patterns confound tissue with chemistry and training size (limits).
- Colon MMRp: MMRp vs MMRd differ in T-cell compartment and in 16 of 43 malignant programmes (Pelka 2021); pooling would average two states; MMRp chosen over MMRd for ≥10 donors and more types with ≥50 cells (CELLxGENE metadata 2026-10-27). Colon results describe MMRp tumours, not CRC in general.
- Lymphoid thoracic only: one homogeneous site, more literature for biological checks; lymph nodes at different sites drain different organs.
- Arm G and chemistry: if the query cannot specify 5′, G may be disadvantaged in the lymph node for technical reasons.
- Partition by donor vs by cell (Spotless): a cell split ignores donor variability and is optimistic; E′ measures that variability.
- Random draw of the held-out donor: typical donor, not the best sampled one. No rotation of the held-out donor (Steyerberg & Harrell 2016) → results hold for the drawn donors (limit, reduced by the held-out group of 4).
- Closed-set rule: unseen test class = dataset shift (Moreno-Torres 2012); Avila Cobos 2020 median RMSE up ≥2-fold in several settings when a type is removed.
- E′ not a ceiling: G averages many donors and can be closer to a new donor than one other donor; P has no such explanation, so P above E′ is examined further.
- Dirichlet α from training: with α = 1 rare types get as many cells as common ones and the anchor reuses rare held-out cells many times; consequence: methods tested less on rare types (Spotless rare-type AUPR not included, BACKLOG).
- Concentration c from Visium HD neighbourhoods: cell types there are coarser than the authors' labels used here (3–6 per tissue); in the lymph node 47% of the nuclei have no type (low UMI per nucleus, similar lymphocyte profiles) and are excluded; the brain sections are mouse, our brain data human. Median, not mean (brain at 10 nuclei, one region c = 110, mean 34).
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
