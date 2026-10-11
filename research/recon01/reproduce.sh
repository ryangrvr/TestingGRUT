#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python research/recon01/audit.py
