"""Pinned low-energy inputs and the canonical three-real-field potential."""
from pathlib import Path
import json,sys,os
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[2]
sys.path.insert(0,str(ROOT.parent/'grut_kinetic_continuation_20261026'/'code'))
import ks_common as K
WORK=Path(os.environ.get('GRUT_ONSET_WORKDIR','/tmp/grut_formation_onset_20261026'));WORK.mkdir(exist_ok=True,parents=True)
MX=K.MX;G=K.GX;LX=K.LX;LP=K.LP;NU=K.V/MX
RC=np.sqrt(LP/G**2)*NU
Q0=K.N_Q0/MX**3;EX=K.RHO_C/MX**4;V0=K.V0/MX**4

def potential(q):
    x,y,z=np.asarray(q);r2=x*x+y*y
    return LP/4*(z*z-NU*NU)**2+G*G/2*z*z*r2+LX/4*r2*r2

def gradient(q):
    x,y,z=np.asarray(q);r2=x*x+y*y;f=G*G*z*z+LX*r2
    return np.array([f*x,f*y,z*(LP*(z*z-NU*NU)+G*G*r2)])

def hessian(q):
    x,y,z=np.asarray(q);r2=x*x+y*y
    return np.array([[G*G*z*z+LX*(r2+2*x*x),2*LX*x*y,2*G*G*x*z],
        [2*LX*x*y,G*G*z*z+LX*(r2+2*y*y),2*G*G*y*z],
        [2*G*G*x*z,2*G*G*y*z,LP*(3*z*z-NU*NU)+G*G*r2]])

def save(name,data):
    def cv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    (ROOT/name).write_text(json.dumps(data,indent=2,default=cv,allow_nan=False)+'\n')

def read(name):return json.loads((ROOT/name).read_text())
