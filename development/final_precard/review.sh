#!/usr/bin/env bash
set -eu
if [ "$#" -ne 2 ]; then
  echo 'Usage: review.sh /absolute/path/GRUT_FINAL_PRECARD_REVIEW_INPUTS.zip /fresh/output/directory' >&2
  exit 2
fi
grut_packet_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$grut_packet_root"
exec "${GRUT_REVIEW_PYTHON:-python3}" -m development.final_precard.run_review --inputs "$1" --output "$2"
