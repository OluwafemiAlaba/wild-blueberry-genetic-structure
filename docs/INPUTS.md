# Required study inputs

Data are distributed separately from this code repository. No public data accession is currently specified. The following paths are relative to the repository root.

## Source files for the full workflow

- `inputs/dose.txt`: tab-separated table containing SNP, CHROM, POS, REF and ALT followed by sample columns. The study source contains 7,407 variant rows and 272 sample columns. Missing calls are NA; observed diploid-style codes are 0, 1 or 2.
- `inputs/amova/angust_geographic.xlsx`: historical V. angustifolium geographic dataset, 135 samples and 6,491 markers.
- `inputs/amova/myrt_regions.xlsx`: 30-sample V. myrtilloides source metadata used for cohort construction.
- `inputs/amova/myrt_without_hope.xlsx`: historical V. myrtilloides geographic dataset, 28 samples and 6,491 markers.
- `inputs/amova/managed_unmanaged.xlsx`: historical 137-sample management dataset. The management script removes Hope34, leaving 136 samples.
- `inputs/manuscript_media/image3.png`: original chromosome-density raster used only to reproduce Figure 2 labeling.

Workbook `Sheet1` stores sample and group in the first two columns, marker headings on row 3 and genotype rows from row 4. Keep source row order and cached numeric values intact. The geographic figure script reads saved `PCoA` coordinates beginning on row 33 and axis percentages in B10/C10. Table export also requires `PhiPT` and `PhiPTP` result sheets.

## Outputs and data ordering

`prepare.py` writes cohort/sample manifests, a variant manifest, `genotypes.npz` and LFMM inputs in `analysis_v2/data/`. LFMM rows follow the analysis sample manifest; columns follow LD-retained variants. LFMM missing calls are encoded as 9. The complete filtered matrix uses NaN for missing calls.

`ancestry.R` writes ancestry proportions, cross-entropy and validation summaries in `analysis_v2/results/`, with fitted/imputed matrices in `analysis_v2/data/`. Ancestry columns have arbitrary component labels, not species probabilities.

`analyze.py` and `clustering.R` write PCA scores, explained-variance values, clustering merges and sample order. Hierarchical-clustering indices follow R's one-based hclust convention. `diagnostics.py` also reads Q matrices saved for individual sNMF runs.

The corrected management matrix, PCoA and AMOVA are written under `outputs/management_Hope34_excluded/`. The pairwise calculation writes `deliverables/Management_differentiation_preview/pairwise_PhiPT_corrected.csv`; that directory name is retained for script compatibility and does not add a separate figure to the manuscript. Table summaries are exported to `outputs/tables/`.

To redraw figures without fitting the analyses, provide the saved output files read by the corresponding figure scripts. The study data and derived output files are not included in this public code repository.
