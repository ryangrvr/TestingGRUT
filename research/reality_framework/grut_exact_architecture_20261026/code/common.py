from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def write(name,data):
    def convert(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,np.ndarray):return x.tolist()
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(data,indent=2,allow_nan=False,default=convert)+'\n')


def read_previous():
    p=ROOT.parent/'grut_heavy_vector_budget_20261025'/'INPUTS.json'
    return json.loads(p.read_text())
