from pathlib import Path
import json,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
WORK=Path(os.environ.get('GRUT_MODE_WORKDIR','/tmp/grut_reheating_modes_20261026'))
WORK.mkdir(exist_ok=True,parents=True)
PARENT=ROOT.parent/'grut_reheating_preparation_20261026'

def write(name,value):
    def native(v):
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError(type(v).__name__)
    (ROOT/name).write_text(json.dumps(value,indent=2,allow_nan=False,default=native)+'\n')
