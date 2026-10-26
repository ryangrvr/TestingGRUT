#!/bin/sh
# kin_source: rebuild every result (about 20 s on one core). Run from this directory.
# Codex files are read only through `git show 9877518:...` / `git show 225028b:...` into the scratchpad.
set -e
cd "$(dirname "$0")"
python3 -B decay_source.py ../results/decay_source_own.json
python3 -B broad_spectrum.py ../results/broad_source.json
python3 -B mediator_source.py ../results/decay_source_own.json ../results/mediator_source.json
python3 -B release_population.py ../results/release_population.json
python3 -B build_interface.py ../results
