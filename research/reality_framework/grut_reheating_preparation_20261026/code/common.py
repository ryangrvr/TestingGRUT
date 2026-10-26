from pathlib import Path
import json,os
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
WORK=Path(os.environ.get('GRUT_REHEATING_WORKDIR','/tmp/grut_reheating_preparation_20261026'))
WORK.mkdir(parents=True,exist_ok=True)


def write(name,value):
    def native(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,np.ndarray):return x.tolist()
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(value,indent=2,allow_nan=False,default=native)+'\n')


def previous():
    return json.loads((ROOT.parent/'grut_heavy_vector_budget_20261025/INPUTS.json').read_text())
