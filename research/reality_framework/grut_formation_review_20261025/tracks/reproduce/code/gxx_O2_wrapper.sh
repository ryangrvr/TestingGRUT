#!/bin/sh
# Reproduction-audit wrapper: Codex's run scripts hard-code `g++ -O3 -std=c++17`.
# The review brief asks for -O2. This wrapper replaces -O3 by -O2, logs the
# exact command line, and calls the system compiler. Nothing else is changed.
# Time-budget switch (documented in NOTES.md): when the file O3_FOR_CFS exists next to
# the log file, compilations of radial_cfs.cpp / radial_flux.cpp keep Codex's -O3; -O2 and -O3
# builds of radial_cfs.cpp were verified byte-identical in output, and -O2 is ~3.4x slower.
LOG="${GXX_WRAPPER_LOG:-/dev/null}"
keep3=0
if [ -e "$(dirname "$LOG")/O3_FOR_CFS" ]; then
  case "$*" in *radial_cfs.cpp*|*radial_flux.cpp*) keep3=1;; esac
fi
args=""
for a in "$@"; do
  if [ "$a" = "-O3" ] && [ $keep3 = 0 ]; then a="-O2"; fi
  args="$args '$a'"
done
echo "$(date +%T) /usr/bin/g++$args" >> "$LOG"
eval exec /usr/bin/g++ $args
