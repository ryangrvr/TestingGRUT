"""Conditional ideal processing with actual finite-Q M1 stationary masses."""
from common import *
import importlib.util
from scipy.integrate import cumulative_trapezoid,simpson

class Seed:
    def __init__(self):
        self.x=PROFILE[:,0];self.p=np.array([OMEGA]);self.y=PROFILE[:,1:].T
        self.j=cumulative_trapezoid(self.x**2*self.y[0]**2,self.x,initial=0)
    def sol(self,r):
        z=[np.interp(r,self.x,v) for v in self.y]
        return np.array(z+[np.interp(r,self.x,self.j)])

def main():
    spec=importlib.util.spec_from_file_location('m1_deformation_source',PARENT/'code/deformation.py')
    D=importlib.util.module_from_spec(spec);spec.loader.exec_module(D)
    sources=json.loads((ROOT/'SOURCE_RESULTS.json').read_text())['rows']
    reference=json.loads((PARENT/'DEFORMATION_RESULTS.json').read_text())['primary']['mass_GeV']
    out=[]
    for tag in ('radial_fine','pulse10','pulse20','pulse40_fine'):
        row=next(r for r in sources if r['tag']==tag)
        free=row['Nplus_out']+row['Nminus_out']+2*row['Nphi_out_net_vacuum_reference']
        q=QREF-free;sol=D.solve(q,.002,Seed(),length=160.,tol=2e-9)
        r=np.linspace(0,160,24001);f,fp,s,sp,j=sol.sol(r);om=float(sol.p[0])
        u=B/4*(s*s-1)**2-C*(s-1)**3+.5*s*s*f*f+A/4*f**4
        e=4*np.pi*simpson(r*r*(.5*(fp*fp+sp*sp)+.5*om*om*f*f+u),x=r)*V/G
        ratio=(e+M*free)/reference
        out.append({'source':tag,'conditional_final_core_charge':q,'free_positive_X':free,
                    'eta':q/QREF,'assumed_final_antiparticles':0.,'actual_stationary_mass_GeV':e,
                    'stationary_E_over_Qm':e/(q*M),'Omega':om,'bound':bool(e<q*M),
                    'cold_mass_over_reference':ratio,'mass_excess_fraction':ratio-1,
                    'within_frozen_one_percent_under_assumptions':bool(abs(ratio-1)<=.01),
                    'max_BVP_residual':float(sol.rms_residuals.max())})
    save('COLD_COMPLETION_RESULTS.json',{'scope':__doc__,'rows':out,
         'actual_capture_and_redshift_demonstrated':False,'actual_final_abundance':None,
         'preparation_or_yield_retuned':False})
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
