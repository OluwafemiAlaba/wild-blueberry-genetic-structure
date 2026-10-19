#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 deliverables/Figure2_revision/relabel_original.py
python3 deliverables/Figure3_revision/create_combined_pcoa.py
python3 deliverables/Figure4_revision/create_figure4.py
python3 deliverables/Figure5_revision/create_management_pcoa.py
