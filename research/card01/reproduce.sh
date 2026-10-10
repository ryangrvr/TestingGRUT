#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python research/card01/audit.py
python research/card01/price_core.py
python research/card01/render_evidence.py
