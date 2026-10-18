#!/usr/bin/env python3
"""Frozen Candidate 1 Q-ball calculation. See results.json for numerical evidence.

Run: python GRUT_Candidate_1/code/qball.py
All dimensionless radii use rho=g*v*r. Scalar conventions: (-+++),
X=exp(-i*omega*t)*f(r)/sqrt(2), F=f/v, S=phi/v.
"""
from __future__ import annotations
import json, platform, sys, time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_bvp, simpson, cumulative_trapezoid
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import diags, bmat, csr_matrix
from scipy.sparse.linalg import eigsh, LinearOperator

OUT=Path(__file__).resolve().parent
CONFIG={
    'a':1.0,'b':0.1,'F0':1.979050307129842,'g':0.1,'v_over_reduced_Mpl':1e-5,
    'Gv2':1e-10/(8*np.pi),
    'boundary_conditions':'Fprime(0)=Sprime(0)=0,F(L)=0,S(L)=1,F(0)=1.979050307129842',
    'bvp_runs':[
        {'name':'base','L':30.,'nodes':301,'tol':1e-5},
        {'name':'grid_refinement','L':30.,'nodes':601,'tol':1e-7},
        {'name':'domain_refinement','L':40.,'nodes':801,'tol':1e-8},
        {'name':'domain_confirmation','L':50.,'nodes':1001,'tol':1e-9},
    ],
    'initial_guesses':[{'width':3.0,'Omega':0.8,'from_preparation_seed':True}],
    'eigen_grids':[{'L':30.,'N':1200},{'L':40.,'N':2000},{'L':40.,'N':3200}],
    'eigen_l':[0,1,2,3],
    'architecture':'Einstein + one complex X + one real phi, exact global U(1)',
    'physical_couplings_frozen_before_results':True,
    'preparation_frozen_before_verification':True,
    'preparation_selection':'Six failed high-central-amplitude BVP starts archived. First fixed-frequency Omega=0.8 bag seed converged, then central F0 was frozen for the verification runs.',
}

def serial(o):
    if isinstance(o,np.ndarray): return o.tolist()
    if isinstance(o,np.generic): return o.item()
    raise TypeError(type(o).__name__)

def log(s):
    print(s,flush=True)
    with (OUT/'results_execution.log').open('a') as f: f.write(s+'\n')

def fun(r,y,p):
    F,Fp,S,Sp=y; O=p[0]
    return np.array([Fp,(S*S+CONFIG['b']*F*F-O*O)*F,Sp,(CONFIG['a']*(S*S-1)+F*F)*S])

def bc(ya,yb,p):
    return np.array([ya[1],ya[3],yb[0],yb[2]-1,ya[0]-CONFIG['F0']])

def seed(r,width):
    F=CONFIG['F0']*np.exp(-(r/width)**2)
    Fp=-2*r/width**2*F
    S=1-.995*np.exp(-(r/(width*.95))**2)
    Sp=2*r/(width*.95)**2*(1-S)
    return np.array([F,Fp,S,Sp])


def main():
    if not (OUT/'results_preparation_seed.npz').exists():
        import runpy
        runpy.run_path(str(OUT/'results_preparation_search.py'),run_name='__main__')
    (OUT/'results_execution.log').write_text('')
    result={'configuration':CONFIG,'versions':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__},'command':'python GRUT_Candidate_1/code/qball.py','attempts':[],'preliminary_failure_archives':['results_initial_failure.json','results_gaussian_failure.json'],'preparation_search':'results_preparation_search.json'}
    (OUT/'results_prescription.json').write_text(json.dumps(result,indent=2,default=serial)+'\n')
    log('Frozen prescription written before numerical solve.')
    sol=None
    for run in CONFIG['bvp_runs']:
        r=np.linspace(0,run['L'],run['nodes'])
        seeds=CONFIG['initial_guesses'] if sol is None else [None]
        accepted=None
        for guess in seeds:
            if guess is None:
                init=sol.sol(r); p=sol.p
                outside=r>sol.x[-1]
                init[:,outside]=np.array([0.,0.,1.,0.])[:,None]
            else:
                seed_data=np.load(OUT/'results_preparation_seed.npz')
                init=np.array([np.interp(r,seed_data['r'],row) for row in seed_data['y']]);p=np.array([guess['Omega']])
            start=time.perf_counter()
            trial=solve_bvp(fun,bc,r,init,p=p,S=np.diag([0.,-2.,0.,-2.]),tol=run['tol'],bc_tol=run['tol'],max_nodes=40000)
            checkr=np.linspace(0,run['L'],10001);checky=trial.sol(checkr)
            attempt={**run,'guess':guess,'status':trial.status,'message':trial.message,'Omega':trial.p[0],'mesh_nodes':len(trial.x),'max_rms_residual':float(np.max(trial.rms_residuals)),'boundary_max_abs':float(np.max(np.abs(bc(trial.y[:,0],trial.y[:,-1],trial.p)))),'min_F':float(checky[0].min()),'min_S':float(checky[2].min()),'F0':float(checky[0,0]),'S0':float(checky[2,0]),'seconds':time.perf_counter()-start}
            result['attempts'].append(attempt);log(json.dumps(attempt,default=serial))
            if trial.success and 0<trial.p[0]<1 and checky[0].min()>-1e-7 and checky[2].min()>0 and checky[2].max()<1.000001:
                accepted=trial;break
        if accepted is None:
            result['error']='No admissible converged BVP solution.'
            (OUT/'results.json').write_text(json.dumps(result,indent=2,default=serial)+'\n')
            return 1
        sol=accepted
        result.setdefault('convergence',[]).append(observables(sol))
    # Remaining computations are separated below for a reproducible refinement pass.
    result['observables']=observables(sol)
    np.savez(OUT/'results_profile.npz',rho=np.linspace(0,50,20001),y=sol.sol(np.linspace(0,50,20001)),Omega=sol.p[0])
    result['spectra']=spectra(sol)
    result['checks']=checks(sol)
    (OUT/'results.json').write_text(json.dumps(result,indent=2,default=serial)+'\n')
    write_csv(sol)
    log('Completed: results.json, results_profile.csv, results_profile.npz')
    return 0


def fields(sol,r):
    F,Fp,S,Sp=sol.sol(r);O=sol.p[0];a=CONFIG['a'];b=CONFIG['b']
    U=a/4*(S*S-1)**2+.5*S*S*F*F+b/4*F**4
    grad=.5*(Fp*Fp+Sp*Sp)
    kinetic=.5*O*O*F*F
    e=grad+kinetic+U
    pr=grad+kinetic-U
    pt=-grad+kinetic-U
    active=e+pr+2*pt
    return {'F':F,'Fp':Fp,'S':S,'Sp':Sp,'U':U,'grad':grad,'kinetic':kinetic,'e':e,'pr':pr,'pt':pt,'active':active}


def observables(sol):
    r=np.linspace(0,sol.x[-1],40001);d=fields(sol,r);O=sol.p[0]
    integrate=lambda x:4*np.pi*simpson(r*r*x,x=r)
    I=integrate(d['F']**2);Q=O*I;E=integrate(d['e'])
    Ecum=4*np.pi*cumulative_trapezoid(r*r*d['e'],r,initial=0)
    r90=np.interp(.9*Ecum[-1],Ecum,r)
    Gv2=CONFIG['Gv2']
    return {'L':sol.x[-1],'Omega':O,'S0':d['S'][0],'Ibar':I,'Qbar':Q,'Ebar':E,'Ebar_over_Qbar':E/Q,'R90_bar':r90,'Egrad_bar':integrate(d['grad']),'Etime_bar':integrate(d['kinetic']),'EU_bar':integrate(d['U']),'integrated_pressure_bar':integrate(d['pr']+2*d['pt']),'active_mass_bar':integrate(d['active']),'Q_at_g_0p1':Q/CONFIG['g']**2,'E_over_v_at_g_0p1':E/CONFIG['g'],'R90_times_v_at_g_0p1':r90/CONFIG['g'],'Schwarzschild_radius_over_R90':2*Gv2*E/r90,'maximum_compactness':np.max(np.divide(2*Gv2*Ecum,r,out=np.zeros_like(r),where=r>0)),'central_energy_density_over_g2v4':d['e'][0],'central_pressure_over_g2v4':d['pr'][0]}


def spectra(sol):
    ans=[];O=float(sol.p[0]);a=CONFIG['a'];b=CONFIG['b']
    for grid in CONFIG['eigen_grids']:
        L=grid['L'];N=grid['N'];h=L/(N+1);r=h*np.arange(1,N+1)
        F,Fp,S,Sp=sol.sol(r)
        # Zero Dirichlet u,w at both ends; centered difference in L2(drho).
        vff=S*S+3*b*F*F-O*O
        vss=a*(3*S*S-1)+F*F
        vfs=2*F*S
        I=observables(sol)['Ibar']
        q=np.concatenate([np.sqrt(4*np.pi)*r*F,np.zeros(N)])
        rank_factor=4*O*O/I
        for ell in CONFIG['eigen_l']:
            cent=ell*(ell+1)/(r*r)
            A=diags([-np.ones(N-1)/h**2,2/h**2+vff+cent,-np.ones(N-1)/h**2],[-1,0,1],format='csr')
            B=diags([-np.ones(N-1)/h**2,2/h**2+vss+cent,-np.ones(N-1)/h**2],[-1,0,1],format='csr')
            C=diags(vfs,0,format='csr')
            H=bmat([[A,C],[C,B]],format='csr')
            if ell==0:
                op=LinearOperator(H.shape,matvec=lambda x: H@x+rank_factor*h*q*np.dot(q,x),dtype=float)
            else: op=H
            start=time.perf_counter()
            ev,vec=eigsh(op,k=5,which='SA',tol=2e-10,maxiter=200000)
            ix=np.argsort(ev);ev=ev[ix];vec=vec[:,ix]
            rec={**grid,'h':h,'ell':ell,'lowest_eigenvalues':ev,'seconds':time.perf_counter()-start,'rank_one_applied':ell==0}
            if ell==1:
                translation=np.concatenate([r*Fp,r*Sp]);translation/=np.linalg.norm(translation)
                rec['lowest_translation_overlap_abs']=abs(float(np.dot(translation,vec[:,0])))
                rec['translation_Rayleigh']=float(translation@(H@translation))
            ans.append(rec);log('spectrum '+json.dumps(rec,default=serial))
    return ans


def checks(sol):
    r=np.linspace(0,sol.x[-1],40001);d=fields(sol,r);O=sol.p[0];Gv2=CONFIG['Gv2']
    integrate=lambda x:4*np.pi*simpson(r*r*x,x=r)
    E=integrate(d['e']);A=integrate(d['active'])
    # Newtonian gauge: ds²=-(1+2 Phi)dt²+(1-2 Psi)dx².
    # Laplacians Phi=4pi G active; Psi=4pi G energy.
    def potential(density):
        m=4*np.pi*cumulative_trapezoid(r*r*density,r,initial=0)
        outer=4*np.pi*(cumulative_trapezoid(r*density,r,initial=0)[-1]-cumulative_trapezoid(r*density,r,initial=0))
        p=-Gv2*(np.divide(m,r,out=np.zeros_like(r),where=r>0)+outer)
        dp=Gv2*np.divide(m,r*r,out=np.zeros_like(r),where=r>0)
        return m,p,dp
    me,Psi,dPsi=potential(d['e']);ma,Phi,dPhi=potential(d['active'])
    # Conservation in flat spacetime: p_r' + 2(p_r-p_t)/rho=0.
    y=sol.sol(r);yp=sol.sol(r,1)
    F,Fp,S,Sp=y;Fpp=yp[1];Spp=yp[3]
    U_F=S*S*F+CONFIG['b']*F**3
    U_S=CONFIG['a']*(S*S-1)*S+F*F*S
    prprime=Fp*Fpp+Sp*Spp+O*O*F*Fp-U_F*Fp-U_S*Sp
    residual=prprime+np.divide(2*(Fp*Fp+Sp*Sp),r,out=np.zeros_like(r),where=r>0)
    scale=max(np.max(np.abs(prprime)),np.max(np.abs(np.divide(2*(Fp*Fp+Sp*Sp),r,out=np.zeros_like(r),where=r>0))))
    grad=integrate(d['grad']);timeE=integrate(d['kinetic']);UE=integrate(d['U'])
    virial=grad+3*(UE-timeE)
    # Phase-direction quadratic operator is positive modulo exact U(1) phase.
    phase=[]
    for grid in CONFIG['eigen_grids']:
        N=grid['N'];h=grid['L']/(N+1);rp=h*np.arange(1,N+1);Fq,_,Sq,_=sol.sol(rp)
        for ell in [0,1,2]:
            diagonal=2/h**2+Sq*Sq+CONFIG['b']*Fq*Fq-O*O+ell*(ell+1)/(rp*rp)
            val=eigh_tridiagonal(diagonal,-np.ones(N-1)/h**2,select='i',select_range=(0,3),check_finite=False)[0]
            phase.append({'ell':ell,'N':N,'h':h,'L':grid['L'],'lowest_eigenvalues':val})
    # Bending at impact parameter b: integral d_z(Phi+Psi) transverse gradient.
    # In dimensionless rho, alpha(B)=2 B int_0^infty [Phi'+Psi'](sqrt(B²+z²))/sqrt(B²+z²) dz.
    from scipy.integrate import quad
    def mass_interp(x,m,total): return np.interp(x,r,m,right=total)
    bends=[]
    for impact in [10.,20.,50.,100.]:
        def integrand(z):
            rr=np.sqrt(impact*impact+z*z)
            return Gv2*(mass_interp(rr,me,E)+mass_interp(rr,ma,A))/rr**3
        alpha=2*impact*quad(integrand,0,np.inf,epsabs=1e-20,epsrel=1e-9,limit=300)[0]
        exterior=4*Gv2*E/impact
        bends.append({'impact_bar':impact,'deflection_rad':alpha,'point_mass_4GME_over_b':exterior,'ratio_to_point_mass':alpha/exterior})
    return {'stress_conservation_max_abs':np.max(np.abs(residual[1:-1])),'stress_conservation_max_relative':np.max(np.abs(residual[1:-1]))/scale,'virial_Tgrad_plus_3_U_minus_Ttime':virial,'virial_relative_to_E':virial/E,'active_mass_over_energy':A/E,'central_Phi':Phi[0],'central_Psi':Psi[0],'central_gravitational_redshift_to_infinity_first_order':-Phi[0],'Phi_at_R90':np.interp(observables(sol)['R90_bar'],r,Phi),'Psi_at_R90':np.interp(observables(sol)['R90_bar'],r,Psi),'phase_sector_spectrum':phase,'light_bending':bends,'asymptotic_F_mass':np.sqrt(1-O*O),'asymptotic_S_mass':np.sqrt(2*CONFIG['a']),'exterior_Phi_at_L':Phi[-1],'exterior_Psi_at_L':Psi[-1],'poisson_active_mass_bar':ma[-1],'poisson_energy_mass_bar':me[-1]}


def write_csv(sol):
    r=np.linspace(0,sol.x[-1],10001);d=fields(sol,r);Gv2=CONFIG['Gv2']
    def potential(x):
        m=4*np.pi*cumulative_trapezoid(r*r*x,r,initial=0)
        cr=cumulative_trapezoid(r*x,r,initial=0)
        return -Gv2*(np.divide(m,r,out=np.zeros_like(r),where=r>0)+4*np.pi*(cr[-1]-cr))
    data=np.column_stack([r,d['F'],d['Fp'],d['S'],d['Sp'],d['e'],d['pr'],d['pt'],d['active'],potential(d['active']),potential(d['e'])])
    np.savetxt(OUT/'results_profile.csv',data,delimiter=',',header='rho,F,Fprime,S,Sprime,energy_over_g2v4,pr_over_g2v4,pt_over_g2v4,active_over_g2v4,Phi,Psi',comments='')

if __name__=='__main__': raise SystemExit(main())
