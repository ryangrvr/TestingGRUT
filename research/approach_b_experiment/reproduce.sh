#!/usr/bin/env bash
set -euo pipefail
packet_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python "$packet_dir/verify_contract.py"
python "$packet_dir/joint_contract.py" \
  --calibration "$packet_dir/../constructive_decision/software_control_calibration.json" \
  --protocol "$packet_dir/../constructive_decision/software_control_protocol.json" \
  --measurement "$packet_dir/software_control_measurement.json" \
  --output "$packet_dir/software_control_joint_contract.json"
