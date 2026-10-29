"""Signed mode-resolved energy flux of the perturbation field, not particles."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def main():
    out={}
    for tag in ('multipoles64','multipoles96','multipoles96_ang12','multipoles96_dr05',
                'multipoles64_l6','multipoles96_l6','multipoles96_t20',
                'multipoles96_double','multipoles96_double_l6','multipoles96_double_dr05','multipoles96_double_t20'):
        z=np.load(ROOT/(tag+'.npz'));a=z['coefficients'];t=z['times'];m=z['modes'];rows=[]
        for ir,r in enumerate((16.,20.,24.)):
            for ch,label in enumerate(('charged_complex','neutral_real')):
                vals=[]
                for ell in range(int(m[:,0].max())+1):
                    select=m[:,0]==ell
                    f=-np.real(np.conj(a[:,ir,1,1,ch,select])*a[:,ir,1,2,ch,select]).sum(axis=1)
                    vals.append(float(np.trapezoid(f,x=t)))
                total=sum(vals)
                rows.append({'radius_mX_inverse':r,'channel':label,'flux_by_l_unscaled':vals,
                             'l2_over_net_flux':vals[2]/total if abs(total)>1e-30 else None,
                             'net_flux_unscaled':total})
        out[tag]={'last_time':float(t[-1]),'rows':rows}
    result={'runs':out,'scope':'Signed time-integrated radial energy flux by l of (disturbed minus undisturbed) classical fields, truncated to listed multipoles. A flux fraction can be negative or exceed one. It omits the interference term with the undisturbed field and is not an outgoing particle population or a relative collision-angle kernel.'}
    (ROOT/'ANGULAR_FLUX.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
