#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 scripts/check_inputs.py
mkdir -p analysis/data analysis/results analysis/Rlib outputs metadata deliverables/Management_differentiation_preview
python3 analysis/prepare.py
Rscript analysis/ancestry.R
python3 analysis/analyze.py
Rscript analysis/clustering.R
python3 scripts/correct_management_hope34.py
python3 scripts/recalculate_pairwise.py
python3 scripts/export_tables.py
sh run_figures.sh
