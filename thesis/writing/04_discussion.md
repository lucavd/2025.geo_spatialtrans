# 4. Discussion

Draft started on 2026-11-25. It holds only the limits already written; the rest follows the results.

## 4.1 Limits

### Cell-type proportions in the colon training data

In the colon dataset, 14 of the 16 training donors also contributed a CD45-enriched aliquot, which accounts for 40% of the training cells (29,308 of 73,062). Because the pseudo-spot compositions and the null baseline are based on the training cell-type proportions, immune cells are over-represented in the colon pseudo-spots relative to the unsorted aliquots. Among the six spot types, tumour epithelial cells are 52% instead of 65%, macrophages 10% instead of 6% and plasma cells 7% instead of 4%. The expression profiles are estimated within cell types, and every arm and the anchor share the same compositions, so this shift affects how closely the spots resemble the tissue but not the comparison between sources.

*[Numbers from `colon_cells.csv.gz` (`scripts/02_split.R`) and `colon_main.rds` (`scripts/03_reference.R`), computed on 2026-11-25; CD45 enrichment in Pelka et al. 2021, doi:10.1016/j.cell.2021.08.003, Methods.]*
