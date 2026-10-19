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
5. **What we need from it:** Set up the ADEMP protocol for our study (`protocol/ADEMP.md`, week 2): aim = do benchmark rankings change with the expression source; DGM = Dirichlet pseudo-spots with NB counts, one per arm; estimand = true cell-type proportions per spot; methods = RCTD, CARD, cell2location, NNLS; performance = RMSE and Jensen–Shannon per method × arm, plus rankings. Report a Monte Carlo SE next to every performance number, and choose the number of replicates from the Monte Carlo SE of the main measure (using the pilot's variability). Set the seed once per study and store the random-number state at the start of each repetition, and **check that the deconvolution packages do not reset the seed internally** (compare the random-number state before and after two runs), otherwise replicates are correlated. If replicates run in parallel, give each process its own random-number stream (e.g. L'Ecuyer-CMRG / rstream) derived from the single seed.
