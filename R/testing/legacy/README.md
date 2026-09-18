# R/testing/legacy — file di baseline pre-S0 (provenienza)

Spostati qui in S0 (2026-09-18) perché **non riproducibili** da `R/testing/full_test.R` così com'era committato (seed fisso 42, n_cells 10 000):

| File | Commit | Config registrata nell'RDS | Prodotto da |
|---|---|---|---|
| `full_test_result.rds` | 8032ec4 (2025-08-06) | n_cells 90 000, pixel_size_um 2.5, griglia 6.5 mm, seed 42 → 180 000 celle, sparsità 86.8 %, UMI medio 9343, 5.21 min | `full_test_visHD.R`, che scriveva sullo **stesso path** di `full_test.R` |
| `full_simulation_data.rds` | 175c6aa (2025-07-31) | n_cells 10 000, pixel_size_um 10, griglia 8 mm, **seed 1200** | `full_test.R` lanciato con `--random-seed` |

Il baseline riproducibile della pipeline attuale è `R/testing/full_test_result_S0.rds` (seed 42, config committata), verificato da `R/testing/test_S0_baseline.R`.
