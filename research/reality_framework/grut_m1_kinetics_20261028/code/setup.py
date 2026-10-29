"""Pinned inputs and serialization. Run with the primary scientific Python."""
from pathlib import Path
import os, sys, json, hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RF = ROOT.parent
SCAT = RF / 'grut_m1_scattering_20261028'
M1 = RF / 'grut_matched_candidate_20261027'
WORK = Path(os.environ.get('GRUT_M1_KINETICS_WORK', str(ROOT / '.work')))
WORK.mkdir(parents=True, exist_ok=True)
os.environ['GRUT_M1_SCATTER_WORK'] = str(WORK / 'inherited')
sys.path.insert(0, str(SCAT / 'code'))
from common import OMEGA, QREF, M, V, G, G2, A, B, C, profile, PROFILE, DYN

MH = np.sqrt(6.)
GH = 2 * G
WIDTH = GH**2 * np.sqrt(1-4/MH**2) / (16*np.pi*MH)
EREF = json.loads((M1/'DEFORMATION_RESULTS.json').read_text())['primary']['mass_GeV']/M

def write_json(name, data):
    def default(x):
        if isinstance(x,np.ndarray): return x.tolist()
        if isinstance(x,np.generic): return x.item()
        if isinstance(x,Path): return str(x)
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(data, indent=2, default=default, allow_nan=False)+'\n')

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
