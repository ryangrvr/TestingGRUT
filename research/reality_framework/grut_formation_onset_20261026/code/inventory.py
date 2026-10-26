"""Extend the earlier necessary thermal-state bound to metastable objects."""
from common import *

def main():
    obj=K.Objects();eta=max(1.,float(max(obj.EQ)));rs=K.RS;rows=[]
    ptemplate=.006571055435898394+.0028047657247771878
    for R in [0.,K.R_VEC['e1e-7_band']+K.R_MED['e1e-7'],K.R_VEC['e1e-7_broad']+K.R_MED['e1e-7'],K.R_VEC['e1e-6_band']+K.R_MED['e1e-6']]:
        e=K.rho_onset(R)/(MX*K.N_Q0)
        for label,p in [('template_inventory',ptemplate),('template_plus_latent_s',ptemplate+1.2787693395),('orbit_inventory',2.5859986321)]:
            old=max(0.,(e-1-rs*p)/(3*(1+2*p)))
            new=max(0.,(e-eta-max(rs,2*eta)*p)/(3*(1+2*p)))
            rows.append(dict(R_total=R,inventory_label=label,p=p,old_minimum_T_over_m=old,extended_minimum_T_over_m=new,
                T_c_over_m=K.Rates().T_c))
    rng=np.random.default_rng(20261026);worst=-np.inf
    for _ in range(10000):
        p=10**rng.uniform(-5,2);nm,ns,nb=p*rng.dirichlet([1,1,1]);nc=rng.uniform(0,1+nm+nb);np_=1+nm+nb-nc
        ep=rng.uniform(.8,eta);em=rng.uniform(.8,eta)
        rest=np_+nm+rs*ns+ep*nc+em*nb;bound=eta+max(rs,2*eta)*p
        worst=max(worst,(rest-bound)/max(bound,1.))
    save('INVENTORY_RESULTS.json',dict(eta=eta,scalar_mass_over_m=rs,rows=rows,random_rest_bound_max_residual=worst,
        inequality='rho/(m nQ) <= eta + max(m_s/m,2 eta)*p + 3 (T/m)(1+2p), eta=max(1,sup E/(m|Q|)).',
        scope='Ideal Bose gases, pressureless objects bounded by supplied eta, no object kinetic/thermal excitations beyond that bound. Covers the tabulated metastable local branch, not arbitrary transient field configurations. Not a formation rule.'))

if __name__=='__main__':main()
