#!/bin/bash
# The workflow sets up its reference data itself in its first step (fetchReferenceData),
# so running this script is no longer needed. It does the same thing by hand: files listed
# in reference_data.json are downloaded unless they are already present in that version.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
python3 "${SCRIPT_DIR}/../bin/fetch_reference_data.py" "${SCRIPT_DIR}/reference_data.json" "${SCRIPT_DIR}"
