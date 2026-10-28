"""Recompute the included cold X6 channels on the CHANGED scalar profile.

The pair term is vacuum Born, not a full in-medium Bogoliubov certificate.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.interpolate import CubicSpline

ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
sys.path.insert(0,str(RF/'grut_cstar_audit_20261024/code'))
sys.path.insert(0,str(RF/'grut_gauge_point_test_20261025/code'))
import charge_leakage as classical
import pair_rate as pairs


class Background:
    def __init__(self):
        self.path=ROOT/'PROFILE_M1.csv';self.hash=hashlib.sha256(self.path.read_bytes()).hexdigest()
        a=np.loadtxt(self.path,delimiter=',',skiprows=1)
        self.f=CubicSpline(a[:,0],a[:,1]);self.s=CubicSpline(a[:,0],a[:,3])
        self.omega=json.loads((ROOT/'DEFORMATION_RESULTS.json').read_text())['primary']['Omega']
    def fields(self,r):
        return self.f(r),self.s(r)
    def potential(self,r,scale=1.):
        f,s=self.fields(r);n=len(np.atleast_1d(r));out=np.zeros((3,3,n))
        out[0,0]=out[1,1]=s*s-1+.2*f*f
        out[2,2]=3*(3*s*s-1)-6*.002/.01*(s-1)+f*f-6
        out[0,1]=out[1,0]=.1*f*f
        out[0,2]=out[2,0]=np.sqrt(2)*s*f
        out[1,2]=out[2,1]=np.sqrt(2)*s*f
        return scale*out


def main():
    bg=Background();coherent=[];born=[]
    for radius,tol,mesh in ((40.,2e-7,1001),(60.,2e-8,1601)):
        row=classical.solve(bg,6,radius,tol,mesh);coherent.append(row)
    for nr,nq in ((6001,96),(12001,192)):
        born.append(pairs.free_pairs(bg,6,nr,nq))
    alpha=json.loads((ROOT/'CONSEQUENCES.json').read_text())['A_GeV_minus2']
    c6=alpha*(2.435e18)**2
    fraction=coherent[-1]['initial_fractional_core_charge_loss_rate_per_second_for_c1']+born[-1]['pair_fractional_core_rate_per_second']
    result={'scope':__doc__,'profile_sha256':bg.hash,'coherent':coherent,'vacuum_Born_pairs':born,
            'matched_alpha6_GeV_minus2':alpha,'Planck_normalized_c6':c6,
            'included_initial_fractional_loss_rate_per_second':fraction*c6*c6,
            'included_initial_rate_times_13_8Gyr':fraction*c6*c6*classical.AGE,
            'complete_quantum_or_hot_population_lifetime':False}
    (ROOT/'COLD_LOSS_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('COLD_INCLUDED_GAMMA_AGE',result['included_initial_rate_times_13_8Gyr'])


if __name__=='__main__':main()
