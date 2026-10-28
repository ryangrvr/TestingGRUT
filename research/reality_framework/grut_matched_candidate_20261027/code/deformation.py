"""One-vacuum scalar-potential deformation and actual stationary-object tests."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_bvp,simpson,cumulative_trapezoid
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
G,LX,LAM,V=.1,.001,.03,.924089637375362


class Stored:
    def __init__(self,q):
        p=RF/f'grut_portal_debris_window_20261023/data/profile_Q{q}_lambda0.03_refined.csv'
        a=np.loadtxt(p,delimiter=',',skiprows=1)
        self.x=a[:,0];self.y=a[:,1:].T
        self.j=cumulative_trapezoid(self.x**2*self.y[0]**2,self.x,initial=0.)
        self.p=np.array([G*G*q/(4*np.pi*self.j[-1])])
    def sol(self,r):
        out=[np.interp(r,self.x,row,left=row[0],right=row[-1]) for row in self.y]
        out.append(np.interp(r,self.x,self.j,left=0.,right=self.j[-1]))
        return np.array(out)


def solve(q,eps,seed,length=120.,tol=2e-8):
    def rhs(r,y,p):
        f,fp,s,sp,j=y
        return np.array([fp,(s*s+LX/G**2*f*f-p[0]**2)*f,
                         sp,LAM/G**2*(s*s-1)*s-3*eps/G**2*(s-1)**2+f*f*s,r*r*f*f])
    def bc(l,r,p):return [l[1],l[3],r[0],r[2]-1,l[4],4*np.pi*p[0]*r[4]-G**2*q]
    x=np.linspace(0,length,1201);guess=seed.sol(np.minimum(x,seed.x[-1]))
    outside=x>seed.x[-1]
    if outside.any():guess[:,outside]=np.array([0,0,1,0,guess[4,~outside][-1]])[:,None]
    sol=solve_bvp(rhs,bc,x,guess,p=seed.p,S=np.diag([0.,-2.,0.,-2.,0.]),tol=tol,bc_tol=tol,max_nodes=40000)
    if not sol.success:raise RuntimeError(sol.message)
    a=sol.sol(np.linspace(0,length,6001))
    if not (0<sol.p[0]<1 and a[0].min()>-1e-7 and a[2].min()>-1e-7):raise RuntimeError('Wrong branch')
    return sol


def analyze(sol,q,eps,tag):
    rho=np.linspace(0,sol.x[-1],24001);f,fp,s,sp,j=sol.sol(rho);om=float(sol.p[0])
    grad=.5*(fp*fp+sp*sp);kin=.5*om*om*f*f
    neutral=LAM/(4*G*G)*(s*s-1)**2-eps/(G*G)*(s-1)**3
    u=neutral+.5*s*s*f*f+LX/(4*G*G)*f**4
    integral=lambda z:float(4*np.pi*simpson(rho*rho*z,x=rho))
    eb=integral(grad+kin+u);mass=eb*V/G
    w=rho*rho*(1-s*s);norm=simpson(w,x=rho)
    def form(qp):return float(simpson(w*np.sinc(qp*rho/(np.pi*G*V)),x=rho)/norm)
    grid=np.linspace(.001,.3,250);ff=[form(x) for x in grid];brackets=[(grid[i],grid[i+1]) for i in range(len(grid)-1) if ff[i]*ff[i+1]<0]
    node=brentq(form,*brackets[0],xtol=1e-13) if brackets else None
    return {'tag':tag,'Q':q,'epsilon':eps,'Omega':om,'mass_GeV':mass,'E_over_Qm':eb/(G*G*q),
            'F0':float(f[0]),'S0':float(s[0]),'radial_domain':float(sol.x[-1]),
            'max_BVP_residual':float(np.max(sol.rms_residuals)),
            'charge_relative_residual':float(om*integral(f*f)/(G*G*q)-1),
            'virial_relative_residual':float((integral(grad)+3*(integral(u)-integral(kin)))/eb),
            'first_form_factor_zero_GeV':node,
            'ideal_node_keV':{name:float(node*node/(2*A*.93149410242)*1e6) for name,A in [('Xe131',131),('Ge74',74),('Ar40',40)]} if node else None,
            'classical_absolute_binding':bool(mass<q*G*V)}


def main():
    rows=[];eps_values=[0.,.001,.0018,.002,.0025,.003]
    seed=Stored(200000);primary=None
    for eps in eps_values:
        sol=solve(200000,eps,seed);seed=sol
        row=analyze(sol,200000,eps,'grid_'+str(eps));rows.append(row)
        if eps==.002:primary=sol
        print('PROFILE',eps,row['E_over_Qm'],row['S0'],row['classical_absolute_binding'],flush=True)
    refined=solve(200000,.002,primary,length=160.,tol=2e-10)
    row=analyze(refined,200000,.002,'primary_refined');rows.append(row)
    r=np.linspace(0,160,3201);z=refined.sol(r)
    np.savetxt(ROOT/'PROFILE_M1.csv',np.column_stack([r,z[:4].T]),delimiter=',',header='rho,F,Fprime,S,Sprime',comments='')
    for q in (199000,201000):
        other=solve(q,.002,refined,length=160.,tol=2e-9)
        rows.append(analyze(other,q,.002,'neighbor_'+str(q)))
    large=solve(1000000,.002,Stored(1000000),length=160.,tol=2e-8)
    rows.append(analyze(large,1000000,.002,'large_charge_control'))
    lo=next(x for x in rows if x['tag']=='neighbor_199000');hi=next(x for x in rows if x['tag']=='neighbor_201000')
    numerical=(hi['mass_GeV']-lo['mass_GeV'])/2000
    expected=row['Omega']*G*V
    result={'scope':'changed classical stationary profiles; not formation, full stability, or full cosmology',
            'potential':'lambda/4*(phi^2-v^2)^2-epsilon*v*(phi-v)^3',
            'primary_epsilon':.002,'single_stationary_vacuum_threshold':LAM*(1-2*np.sqrt(2)/3),
            'new_coupling_is_supplied':True,'primary':row,'rows':rows,
            'first_law':{'neighbor_dE_dQ_GeV':numerical,'Omega_mX_GeV':expected,'relative_difference':numerical/expected-1,
                         'neighbor_dOmega_dQ':(hi['Omega']-lo['Omega'])/2000},
            'domain_wall_claim':'one stationary vacuum removes the two-vacuum phi domain-wall mechanism; no all-defect or nonequilibrium-interface theorem'}
    (ROOT/'DEFORMATION_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
