# Esempio 1: dataset "facile" (ben separato, con marker forti)
simulate_spatial_transcriptomics(
  image_path = here::here("images/colon.png"),
  output_path = "data/simulated_easy_correlation.rds",
  output_plot = "results/simulated_easy_distribution.png",
  difficulty_level = "easy",
  n_cells = 20000,
  n_genes = 100
)


# Esempio: genera grafici diagnostici per un dataset simulato
generate_diagnostic_plots(
  input_rds_path = "data/simulated_easy_correlation.rds",
  output_dir = "images/simulation_plots",
  base_name = "sim_analysis"  # opzionale, altrimenti usa il nome del file
)


analyze_and_compare_clusters(rds_path = "data/simulated_easy_correlation.rds")

run_simulation_testing("data/simulated_easy_correlation.rds", "results/simulation_tests/")
