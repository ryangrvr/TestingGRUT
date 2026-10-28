"""Limited nonresonant vacuum 2-to-2 kernels; not a full plasma collision operator.

External phi states use the narrow-width approximation. The resonant X Xbar
elastic/inverse-decay overlap and in-medium/multiparticle processes stay open.
"""
from common import *
from scipy.integrate import quad

def momentum(s,m):return np.sqrt(max(0,s/4-m*m))

def matrix_pair(s,costheta,eps):
    px=momentum(s,1.);ph=momentum(s,np.sqrt(6))
    t=1+6-s/2+2*px*ph*costheta;u=2+12-s-t
    lam=.2;cubic=6*(.03-eps)/G
    return .02+lam*cubic/(s-6)+lam*lam*(1/(t-1)+1/(u-1))

def pair_cross_sections(s,eps):
    px=momentum(s,1.);ph=momentum(s,np.sqrt(6))
    if px==0 or ph==0:return 0.,0.
    angular=quad(lambda z:matrix_pair(s,z,eps)**2,-1,1,epsabs=1e-15,epsrel=2e-12)[0]
    # Identical final phi pair gets 1/2, distinguishable X Xbar does not.
    return angular/(64*np.pi*s)*ph/px,angular/(32*np.pi*s)*px/ph

def neutral_elastic(s,eps):
    p=momentum(s,np.sqrt(6));g3=6*(.03-eps)/G
    def integrand(z):
        t=-2*p*p*(1-z);u=-2*p*p*(1+z)
        amp=.18+g3*g3*(1/(s-6)+1/(t-6)+1/(u-6))
        return amp*amp
    return quad(integrand,-1,1,epsabs=1e-15,epsrel=2e-12)[0]/(64*np.pi*s)

def main():
    rows=[]
    for ecm in (4.95,5.,6.,8.,12.,20.):
        s=ecm*ecm;new,reverse=pair_cross_sections(s,.002);old,_=pair_cross_sections(s,0.)
        hn=neutral_elastic(s,.002);ho=neutral_elastic(s,0.)
        px=momentum(s,1.);ph=momentum(s,np.sqrt(6))
        rows.append({'sqrt_s_over_mX':ecm,'sigma_XbarX_to_phiphi_mXminus2':new,
                     'sigma_phiphi_to_XbarX_mXminus2':reverse,'pair_kernel_new_over_old':new/old,
                     'sigma_phiphi_elastic_mXminus2':hn,'neutral_elastic_new_over_old':hn/ho,
                     'detailed_balance_relative_error':2*px*px*new/(ph*ph*reverse)-1})
    out={'scope':__doc__,'rows':rows,'updated_full_2_to_2_operator':False,'updated_4_to_2_operator':False,
         'open':'Bose/in-medium corrections, real-intermediate-state subtraction, scalar reservoir chemistry and high-multiplicity kernels'}
    save('VACUUM_KERNELS.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__':main()
