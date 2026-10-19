#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 scripts/check_inputs.py
mkdir -p analysis_v2/data analysis_v2/results analysis_v2/Rlib outputs metadata deliverables/Management_differentiation_preview
python3 analysis_v2/prepare.py
Rscript analysis_v2/ancestry.R
python3 analysis_v2/analyze.py
python3 analysis_v2/diagnostics.py
Rscript analysis_v2/clustering.R
python3 analysis_v2/extra_checks.py
python3 scripts/correct_management_hope34.py
python3 scripts/recalculate_pairwise.py
python3 scripts/export_tables.py
sh run_figures.sh
