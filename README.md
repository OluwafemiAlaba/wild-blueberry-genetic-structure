# Wild-blueberry genetic diversity and population structure

Analysis scripts for the manuscript *SNP-based assessment of geographic and management-associated genetic structure in wild blueberries in Maine*.

**Manuscript status:** Under review.

## Requirements

Python 3.12 with the packages in `requirements.txt`; R with **LEA** and **cluster**. The recorded R environment used R 4.6.0 and LEA 3.24.0.

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Data and execution

This repository contains code only. Study data are distributed separately; a public archive accession has not yet been assigned. Required files, formats and output paths are described in [docs/INPUTS.md](docs/INPUTS.md).

Run commands from the repository root. To run the complete workflow and check its outputs:

```sh
sh run_analysis.sh
python3 scripts/validate_package.py
```

The workflow overwrites derived outputs. Individual steps are listed below in execution order.

## 1. Genotype filtering

Construct the sample cohort, filter biallelic SNPs and prune linkage disequilibrium. The study retains 129 samples and 4,640 SNPs after filtering and LD pruning.

```sh
python3 scripts/check_inputs.py
mkdir -p analysis_v2/data analysis_v2/results analysis_v2/Rlib outputs metadata
python3 analysis_v2/prepare.py
```

## 2. Ancestry estimation and imputation

Fit LEA sNMF models with five runs per K: K = 1–6 for all retained samples and K = 1–5 for V. angustifolium. Evaluate cross-entropy, assess imputation using held-out calls and impute missing genotypes while retaining observed calls.

```sh
Rscript analysis_v2/ancestry.R
```

## 3. Ordination and clustering

Calculate standardized PCA and Ward.D2 hierarchical clustering. Additional checks examine missingness, technical replicates, observed-call sensitivity and within-species structure.

```sh
python3 analysis_v2/analyze.py
python3 analysis_v2/diagnostics.py
Rscript analysis_v2/clustering.R
python3 analysis_v2/extra_checks.py
```

## 4. Management-associated differentiation

Exclude Hope34, identified as V. myrtilloides, from the V. angustifolium management dataset. Recalculate PCoA, AMOVA and pairwise PhiPT, then export table summaries.

```sh
python3 scripts/correct_management_hope34.py
python3 scripts/recalculate_pairwise.py
python3 scripts/export_tables.py
```

The geographic and management analyses use the historical 6,491-marker datasets. These differ from the filtered dataset used in steps 1–3.

## 5. Figures

Generate Figures 2–5 from the required inputs and analysis outputs:

```sh
sh run_figures.sh
```

## Using other datasets

Update input paths, workbook layouts and sample identifiers to match your data. In `analysis_v2/prepare.py`, review cohort assignments, missingness and allele-frequency thresholds, and LD settings. In `analysis_v2/ancestry.R`, review K ranges and fitting settings. Update study-specific sample counts in validation scripts and group labels in figure scripts. Keep sample and marker order consistent throughout. Random seeds are specified in the analysis scripts.

The study uses diploid-style 0/1/2 coding; these values are not measured tetraploid allele dosages. Input conventions and historical processing limits are described in [docs/INPUTS.md](docs/INPUTS.md).

## Software

Ancestry estimation uses LEA sNMF. The upstream genotyping pipeline is [GBSapp](https://github.com/bodeolukolu/GBSapp). Python dependencies are listed in `requirements.txt`.

No reuse license has been assigned yet.

## Contact

Oluwafemi A. Alaba: [oluwafea@ucr.edu](mailto:oluwafea@ucr.edu)
