#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python research/recon02/audit_routing.py
