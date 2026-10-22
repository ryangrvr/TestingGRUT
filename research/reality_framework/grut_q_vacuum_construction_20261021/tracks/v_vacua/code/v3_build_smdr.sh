#!/bin/sh
# Download and build SMDR v1.3 (Martin-Robertson, arXiv:1907.02500) in the scratchpad, then compile our driver.
# Nothing of SMDR is placed in the repository.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
DL=${SMDR_DL:-/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/v_vacua_dl/smdr}
mkdir -p "$DL"
cd "$DL"
[ -f smdr-1.3.tar.gz ] || curl -sSL -o smdr-1.3.tar.gz https://www.niu.edu/spmartin/smdr/smdr-1.3.tar.gz
[ -d SMDR-1.3 ] || tar xzf smdr-1.3.tar.gz
cd SMDR-1.3
[ -f libsmdr.a ] || make > ../build.log 2>&1
gcc -O2 -fcommon -I. -o smdr_driver "$HERE/smdr_driver.c" -L. -lsmdr -l3vil -ltsil -lm
echo "built $DL/SMDR-1.3/smdr_driver"
