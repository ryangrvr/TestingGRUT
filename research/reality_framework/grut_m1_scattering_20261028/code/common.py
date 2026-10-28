"""Fixed M1 point, profile and outputs for the antiparticle scattering gate."""
from pathlib import Path
import os,json,hashlib,subprocess,time
import numpy as np
from scipy.interpolate import CubicHermiteSpline

ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
PARENT=RF/'grut_matched_candidate_20261027';DYN=RF/'grut_m1_dynamics_20261027'
WORK=Path(os.environ.get('GRUT_M1_SCATTER_WORK','/workspace/scratch/0c6d6cf183c6/m1_scattering_runs'))
WORK.mkdir(parents=True,exist_ok=True)
OLD_WORK=Path(os.environ.get('GRUT_M1_DYNAMICS_WORK','/workspace/scratch/0c6d6cf183c6/m1_dynamics_runs'))
POINT=json.loads((PARENT/'CANDIDATE_M1.json').read_text())
PROFILE=np.loadtxt(PARENT/'PROFILE_M1.csv',delimiter=',',skiprows=1)
OMEGA=json.loads((PARENT/'DEFORMATION_RESULTS.json').read_text())['primary']['Omega']
G=.1;G2=.01;A=.1;B=3.;C=.2;QREF=200000.;M=POINT['light_renormalized_targets']['mX_GeV']
V=POINT['light_renormalized_targets']['vPhi_GeV'];MPHI=np.sqrt(6)*M
FS=CubicHermiteSpline(PROFILE[:,0],PROFILE[:,1],PROFILE[:,2])
SS=CubicHermiteSpline(PROFILE[:,0],PROFILE[:,3],PROFILE[:,4])

def profile(r):
    r=np.asarray(r);f=FS(np.minimum(r,160));s=SS(np.minimum(r,160))
    return np.where(r>160,0,f),np.where(r>160,1,s)

def save(name,data):
    def default(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(data,indent=2,default=default,allow_nan=False)+'\n')

