#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 outputs/Figure2_revision/relabel_original.py
python3 outputs/Figure3_revision/create_combined_pcoa.py
python3 outputs/Figure4_revision/create_figure4.py
python3 outputs/Figure5_revision/create_management_pcoa.py
