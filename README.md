# SNP-based assessment of geographic and management-associated genetic structure in wild blueberries in Maine

R and Python scripts for genotype filtering, population structure analyses, ordination, hierarchical clustering, AMOVA summaries, and manuscript figures.

This repository contains **code and documentation only**. Genotype data, sample metadata, saved numerical results, manuscript text and photographs are not distributed here. A public data archive accession has not yet been assigned. Running the study analyses requires the source files described in [docs/INPUTS.md](docs/INPUTS.md).

## Installation

Use Python 3.12 and the dependencies in `requirements.txt`:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
```

The ancestry analysis requires R and the Bioconductor package **LEA**; hierarchical clustering additionally uses **cluster**. The recorded analysis environment was R 4.6.0 and LEA 3.24.0. Install a compatible Bioconductor release and record any version changes.

## Running the analyses

Place the source files in the directory structure described in `docs/INPUTS.md`, then run:

```sh
python3 scripts/check_inputs.py
sh run_analysis.sh
python3 scripts/validate_package.py
```

The main workflow performs genotype QC and LD pruning, sNMF ancestry estimation and imputation, standardized PCA, Ward.D2 clustering, diagnostic analyses, the management-cohort PCoA and AMOVA, pairwise differentiation, table exports and figure generation. Running it replaces derived outputs in the working directory.

To redraw figures from existing results:

```sh
sh run_figures.sh
```

## Script organization

- `analysis_v2/prepare.py`: genotype QC, cohort construction and LD pruning.
- `analysis_v2/ancestry.R`: repeated sNMF fits, cross-entropy evaluation, imputation validation and final imputation.
- `analysis_v2/analyze.py`: standardized PCA, kernel-based tests, sensitivity analyses and technical-replicate comparisons.
- `analysis_v2/clustering.R`: Euclidean-distance hierarchical clustering using Ward.D2.
- `analysis_v2/diagnostics.py` and `extra_checks.py`: diagnostic and sensitivity analyses. These are not all main manuscript results.
- `scripts/correct_management_hope34.py`: removal of the author-identified V. myrtilloides sample from the V. angustifolium management analysis; recalculation of PCoA and AMOVA.
- `scripts/recalculate_pairwise.py`: pairwise PhiPT for the corrected management cohort.
- `scripts/export_tables.py`: machine-readable main and supplementary table summaries.
- `scripts/validate_package.py`: dataset dimensions, sample alignment and observed-call preservation checks.
- `deliverables/Figure2_revision/relabel_original.py`: labels an existing SNP-density raster without changing the plotted data.
- `deliverables/Figure3_revision/create_combined_pcoa.py`: two species-specific geographic PCoAs from saved GenAlEx coordinates.
- `deliverables/Figure4_revision/create_figure4.py`: cross-entropy, clustering, ancestry proportions and PCA.
- `deliverables/Figure5_revision/create_management_pcoa.py`: corrected management PCoA.

## Interpretation and reproducibility

The historical 6,491-marker analyses and the newer 129-sample, 4,640-marker analyses use different datasets. Ancestry estimates in the newer workflow are obtained with **LEA sNMF**, not the STRUCTURE or ADMIXTURE programs. Seeds and numerical settings are specified in the scripts.

The supplied 0/1/2 encoding is diploid-style coding and should not be interpreted as measured tetraploid allele dosage. Upstream recoding and missing-call handling for the historical matrices are not reconstructed by these scripts. The newer sNMF imputation is documented separately and preserves observed calls.

## Software attribution

The scripts depend on LEA, R and the Python libraries listed in `requirements.txt`. The upstream sequencing/genotyping pipeline is described in [GBSapp](https://github.com/bodeolukolu/GBSapp); its software is not redistributed here. The exact upstream commands and reference assembly version must be obtained from the study's sequencing records.

No reuse license has been assigned to this repository yet.
