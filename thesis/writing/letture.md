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
