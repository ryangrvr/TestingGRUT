"""Read-only access to the inherited kinetic packet; new output directory."""
from pathlib import Path
import os, sys, json, hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RF = ROOT.parent
KIN = RF / 'grut_m1_kinetics_20261028'
WORK = ROOT / '.work'
WORK.mkdir(exist_ok=True)
os.environ['GRUT_M1_KINETICS_WORK'] = str(WORK / 'inherited')
sys.path.insert(0, str(KIN / 'code'))
import setup as fixed

def save(name, obj):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(obj, indent=2, default=convert, allow_nan=False)+'\n')

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
