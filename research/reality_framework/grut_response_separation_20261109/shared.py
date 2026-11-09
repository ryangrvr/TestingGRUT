import sys,json,importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent
D1=HERE.parent/'grut_source_dictionary_20261109';B2=HERE.parent/'grut_detector_reconstruction_20261108';N22=HERE.parent/'grut_n22_refinement_20261109'
sys.path.insert(1,str(B2));sys.path.insert(2,str(N22))
from model import reconstruction
from lattice import schwinger
from resource_limit import bounded_run
spec=importlib.util.spec_from_file_location('d1_reference',D1/'reference.py');d1_reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(d1_reference)
def card():return json.loads((HERE/'protocol.json').read_text())
def write(path,value):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def encode(matrix):return [[[float(v.real),float(v.imag)] for v in row] for row in matrix]
