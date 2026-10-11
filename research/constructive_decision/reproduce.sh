#!/usr/bin/env bash
set -euo pipefail
packet_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python "$packet_dir/verify.py"
python "$packet_dir/forecast.py" \
  --calibration "$packet_dir/software_control_calibration.json" \
  --protocol "$packet_dir/software_control_protocol.json" \
  --output "$packet_dir/software_control_forecast.json"
