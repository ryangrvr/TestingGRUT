#!/usr/bin/env python3
"""Conditional high-ell subgap comparison, evaluated on the saved profile."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
REF=HERE.parent/'grut_candidate_1_20261018/code/results_profile.npz'


def main():
    start=time.perf_counter()
    protocol=(HERE/'angular_bound_protocol.json').read_bytes()
    config=json.loads(protocol)
    data=np.load(REF);rho=data['rho'];spline=CubicSpline(rho,data['y'],axis=1)
    def integrand(r):
        F,_,S,_=spline(r)
        a=S*S+.3*F*F-1;b=S*S+.1*F*F-1
        d=3*(S*S-1)+F*F;m=2*F*S
        least=np.minimum(b,(a+d-np.sqrt((a-d)**2+4*m*m))/2)
        return np.asarray(r)**2*np.maximum(0,-least)
    grids=[];best=None
    for n in config['numerics']['grids_on_saved_profile_domain']:
        r=np.linspace(float(rho[0]),float(rho[-1]),n);values=integrand(r);j=int(np.argmax(values))
        grids.append({'points':n,'sampled_B':float(values[j]),'rho_at_sampled_max':float(r[j])})
        bracket=(r[max(0,j-1)],r[min(n-1,j+1)])
        opt=minimize_scalar(lambda t:-float(integrand(t)),bounds=bracket,method='bounded',options={'xatol':1e-12})
        best={'optimized_B':float(-opt.fun),'rho_at_max':float(opt.x),'success':bool(opt.success),'bracket':list(map(float,bracket))}
    Bsafe=int(np.ceil(best['optimized_B']))+2
    ell=2
    while ell*(ell+1)<=Bsafe:ell+=1
    tails=[]
    r=np.linspace(float(rho[0]),float(rho[-1]),16001)
    for cutoff in [10,15,20,25]:
        mask=r>=cutoff
        if mask.any():tails.append({'rho_min':cutoff,'maximum_sampled_B_integrand':float(np.max(integrand(r[mask])))})
    result={'protocol_sha256':hashlib.sha256(protocol).hexdigest(),'profile_sha256':hashlib.sha256(REF.read_bytes()).hexdigest(),
            'grid_results':grids,'optimized':best,'saved_domain':[float(rho[0]),float(rho[-1])],
            'numerical_margin_comparison_value':Bsafe,'conditional_first_excluded_ell':ell,
            'tail_diagnostics':tails,'endpoint_fields_F_S':[float(spline(rho[-1])[0]),float(spline(rho[-1])[2])],
            'half_line_bound_certified':False,'scope':'Conditional angular comparison and numerical diagnostic; not a certified entire continuum spectrum or determinant.',
            'elapsed_seconds':time.perf_counter()-start}
    (HERE/'angular_bound_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
