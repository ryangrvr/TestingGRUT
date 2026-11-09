import sys,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
B2=HERE.parent/'grut_detector_reconstruction_20261108'
N22=HERE.parent/'grut_n22_refinement_20261109'
sys.path.insert(1,str(B2));sys.path.insert(2,str(N22))
from model import reconstruction
from lattice import schwinger
from resource_limit import bounded_run
def card():return json.loads((HERE/'protocol.json').read_text())
def write(path,value):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
