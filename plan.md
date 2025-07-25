# Spatial Transcriptomics Simulation Improvement Plan

## Notes
- Main simulation file is R/run_full_size_optimized.R, which includes a biological validation step.
- The goal is to generate spatial transcriptomics data that is both realistic and suitable for simulations.
- Current issue: k-means clustering produces unrealistic round clusters.
- User is seeking improved clustering approaches (see CLUSTER_IMPROV.md) and project context (see CLAUDE.md).
- User wants to focus on improving DBSCAN+Graph pipeline, with minimal and simple code changes.
- Validation should include a check for improved cluster morphology and output a summary image similar to the provided example.
- DBSCAN+Graph pipeline is already in use in the main script.
- After clustering, NA values in intensity_cluster are now filtered out to avoid downstream errors.
- Quick script for cluster statistics after simulation has been added (R/quick_cluster_stats.R).
- User approved adding cluster morphology validation and plot to the biological validation report.
- The validation script should not source or re-run the simulation pipeline; this needs correction.
- Error: validate_cluster_morphology function was missing; needs to be implemented as a local function in the validation script.
- Patch applied: cluster morphology validation function is now implemented locally and the report checks n_noise directly (no $noise access). Next, verify that sparse matrix warnings are resolved.
- Validation for large clusters now samples points to avoid memory warnings; validation is efficient for large datasets.
- Warning: sparse-to-dense coercion warning still present; further optimization needed.
- User requested to extend sparse-to-dense avoidance to all validation functions, not just validate_spatial_coherence.
- plot_distances_by_cluster has now been optimized to avoid unnecessary sparse-to-dense coercion.
- Persistent memory error: fallback to dense path in validate_expression_ranges (apply(..., 1, max)) is no longer needed.
- User requested an additional cluster morphology image: silhouette/contour of cluster shapes for visual inspection.
- User requested geometric metrics (e.g., area, eccentricity) for each cluster in addition to the silhouette plot.
- Validation must always save cluster metric CSV and silhouette/morphology PNG files, even if results are already present in memory or on disk.
- Silhouette plot output currently shows excessive overlap, indicating poor cluster separation or visualization; next steps should address this issue.
- User requested a full theoretical review of the analysis pipeline, based on README.md and CLAUDE.md.
- New direction: generate a synthetic tissue image (with controllable complexity) and use graph-based clustering on expression profiles to obtain realistic, irregular clusters for simulation benchmarking. This approach balances simplicity and biological plausibility, as discussed and agreed.
- Synthetic tissue image generation (Voronoi/patches) function was debugged and now produces correct, non-uniform output. Function tested and validated with console output and PNG.
- Synthetic tissue image integration into the main simulation pipeline implemented and tested; normalization and thresholding fixed.

## Task List
- [x] Analyze R/run_full_size_optimized.R, focusing on the biological validation and clustering steps.
- [x] Review CLUSTER_IMPROV.md for proposed clustering improvements.
- [x] Review CLAUDE.md to understand project goals and requirements.
- [x] Identify and suggest improvements for more realistic spatial clustering and validation.
- [x] Modify clustering step to use/optimize DBSCAN+Graph pipeline with minimal code changes.
- [x] Add validation to produce a summary cluster morphology image (as in the example).
- [x] Run the new pipeline and check if cluster shapes are improved and validation image is generated.
- [x] Review and update R/biological_validation_report.R for correctness and alignment with new pipeline.
  - [x] Add cluster morphology validation and plot to the report (implement as a local function if missing).
  - [x] Remove any code that sources or re-runs the simulation pipeline during validation.
- [x] Test the validation and ensure no sparse-to-dense conversion warnings remain.
- [ ] Investigate and eliminate remaining sparse-to-dense conversion warnings in validation.
  - [x] Audit all validation functions for unnecessary sparse-to-dense coercion and optimize as needed in plot_distances_by_cluster.
  - [x] Audit remaining validation functions for unnecessary sparse-to-dense coercion and optimize as needed, especially validate_expression_ranges (remove apply(..., 1, max) fallback for dense matrices).
  - [x] Add a silhouette/contour plot of cluster shapes as an additional validation image.
  - [x] Compute and save geometric metrics (area, eccentricity, etc.) for each cluster as part of validation.
  - [ ] Improve cluster separation and/or visualization for meaningful silhouette output.
  - [x] Review and summarize the theoretical analysis pipeline from README.md and CLAUDE.md.
  - [x] Design and test synthetic tissue image generation function (Voronoi/patches) and validate output.
  - [x] Integrate synthetic image option (use_synthetic_image) into main simulation pipeline and test end-to-end flow.
  - [ ] Assign initial regions, simulate expression, perform graph-based clustering on expression to obtain realistic clusters, validate morphology/metrics as before.

## Current Goal
- Assign initial regions, simulate expression, perform graph-based clustering on expression to obtain realistic clusters, validate morphology/metrics as before.