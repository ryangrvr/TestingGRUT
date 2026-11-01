"""fx_chain replacement of Codex's bootstrap.py: read the kin_* interfaces from the review tree in place (they are
unchanged since e71730c: checked by m1_kin.py against the pinned SOURCES.json hashes) and cache in the scratchpad."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[2] / 'results'
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
WORK = Path(os.environ.get('GRUT_KINETIC_WORKDIR', SP + '/wt2/fx_chain/kwork'))
WORK.mkdir(parents=True, exist_ok=True)
INPUT_ROOT = Path(os.environ.get('GRUT_KINETIC_INPUT_ROOT',
                                 '/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025'))
