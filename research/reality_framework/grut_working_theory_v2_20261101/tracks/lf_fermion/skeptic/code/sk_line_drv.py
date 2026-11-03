"""Skeptic driver: runs lf_fermion's own lfline.main (read-only import) with the output redirected to the skeptic's
scratch (lf_fermion's recorded runs are only read). Used for the N1 fermion-grid check (Nf 8192 -> 16384) that
lf_fermion did only for a chi line. Claude, lf_fermion skeptic (issue #13 context).
Usage: python3 -B sk_line_drv.py <lfline arguments>"""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'code'))
import lfcommon as LC
LC.SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/lf_fermion_skeptic')
import lfline
sys.argv = ['lfline.py'] + sys.argv[1:]
sys.exit(lfline.main())
