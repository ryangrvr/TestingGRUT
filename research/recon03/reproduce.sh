#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/recon03/audit.py
