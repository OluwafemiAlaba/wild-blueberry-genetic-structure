"""Check required source files before starting the analysis."""
from pathlib import Path
import sys
root = Path(__file__).resolve().parents[1]
required = ["inputs/dose.txt", "inputs/amova/angust_geographic.xlsx",
            "inputs/amova/myrt_regions.xlsx", "inputs/amova/myrt_without_hope.xlsx",
            "inputs/amova/managed_unmanaged.xlsx", "inputs/manuscript_media/image3.png"]
missing = [p for p in required if not (root / p).is_file()]
if missing:
    print("Required study inputs are not included in this code-only repository:", file=sys.stderr)
    for path in missing:
        print("  " + path, file=sys.stderr)
    print("See docs/INPUTS.md for file formats and directory placement.", file=sys.stderr)
    sys.exit(1)
print("All required source files are present.")
