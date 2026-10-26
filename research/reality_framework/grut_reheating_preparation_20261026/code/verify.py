"""Independent controls and mathematical checks for the actual new runs."""
import json,math
import numpy as np
from scipy.integrate import solve_ivp
from common import ROOT,WORK,write,previous
from reheating import build,run,rhs,geometric,potential_force,XI,ETA2,COLS


def main():
    data=json.loads((ROOT/'REHEATING_RESULTS.json').read_text());rows={x['tag']:x for x in data['rows']}
    central=rows['b2_g0.1'];exe=build()
    short,_=run(exe,'short_control',2.,.1,end=10.,samples=1000,tol=2e-12)
    sol=solve_ivp(lambda t,v:rhs(t,v,2.,.1,XI,ETA2),(0,10),short['initial_state'],method='Radau',rtol=2e-10,atol=2e-13)
    assert sol.success,sol.message
    c=np.array(short['final_state']);s=sol.y[:,-1]
    independent=float(np.max(abs(c-s)/(1+abs(c))))
    long,arrlong=run(exe,'long_control',2.,.1,end=500.,samples=10000,tol=2e-11)
    write('INDEPENDENT_REHEATING.json',dict(short_C=short,Radau_endpoint=s.tolist(),Radau_nfev=sol.nfev,Radau_state_error=independent,long_control=long))

    # Solve trace, Raychaudhuri and F'' simultaneously, independently of
    # the trace-elimination expression used by the numerical equations.
    rng=np.random.default_rng(262026);errors=[];wards=[]
    for _ in range(50):
        v=np.zeros(19);v[0]=rng.uniform(.1,2);v[1]=rng.normal();v[14]=rng.uniform(0,3)
        z=np.array([rng.normal()+1j*rng.normal(),rng.normal()+1j*rng.normal(),28.3+1j*rng.normal()])
        vel=rng.normal(size=3)+1j*rng.normal(size=3)
        v[2:8:2]=z.real;v[3:8:2]=z.imag;v[8:14:2]=vel.real;v[9:14:2]=vel.imag
        eta=10**rng.uniform(-9,-5);pot,F=potential_force(z);K=np.vdot(vel,vel).real
        f=1-2*XI*eta*abs(z[0])**2;fd=-4*XI*eta*np.real(z[0].conjugate()*vel[0])
        rho=.5*v[1]**2+2*v[0]**2+eta*(K+pot)+v[14]
        v[16]=(-fd+math.sqrt(fd*fd+4*f*rho/3))/(2*f)
        z,vel,pot,F,f,fd,R,hd=geometric(v,2.,XI,eta)
        S=v[1]**2+2*eta*K+4*v[14]/3
        base=-4*XI*eta*(abs(vel[0])**2-3*v[16]*np.real(z[0].conjugate()*vel[0])-np.real(z[0].conjugate()*F[0]))
        mat=np.array([[1.,-6.,0.],[0.,2*f,1.],[-4*XI*XI*eta*abs(z[0])**2,0.,1.]])
        expected=np.linalg.solve(mat,np.array([12*v[16]**2,-S+v[16]*fd,base]))
        errors.append(max(abs(expected[0]-R)/(1+abs(R)),abs(expected[1]-hd)/(1+abs(hd))))
        ward=np.dot([1,3,6],np.imag(z.conjugate()*F))
        wards.append(abs(ward)/max(1.,sum(abs(z.conjugate()*F))))
    # Analytic early matter-origin growing mode: X=t^(2/3), H=2/(3t),
    # xi R=-10/(9t^2). It is an origin test, not the charged condensate.
    def growth(t,y):return [y[1],-2*y[1]/t+10*y[0]/(9*t*t)]
    grow=solve_ivp(growth,(1,1000),[1.,2/3],method='DOP853',rtol=1e-11,atol=1e-13)
    growth_error=abs(grow.y[0,-1]/100-1)
    old=json.loads((ROOT.parent/'grut_heavy_vector_budget_20261025/BACKGROUND_RESULTS.json').read_text())['late']
    rd=rows['RD_reference']['reference_window_30_40']
    rd_errors=dict(NX=abs(rd['NX']/old['N_X']-1),EX=abs(rd['EX']/old['E_X']-1),epsilon=abs(rd['epsilon']/old['epsilon']-1),
                   NY=abs((rd['NYH']+rd['NYL'])/(old['N_YH']+old['N_YL'])-1))
    tight=rows['central_tight'];tolerr=max(abs(central[k]/tight[k]-1) for k in ('late_NX','late_EX','late_epsilon'))
    phase_slope=(rows['phase_plus']['late_NX']-rows['phase_minus']['late_NX'])/(2e-8*abs(central['late_NX']))
    checks={k:max(x['checks'][k] for x in rows.values()) for k in ('max_Friedmann','max_Gauss','max_trace','max_total_work','max_hidden_work')}
    phase_response=max(abs(rows[t]['late_NX']/central['late_NX']-1) for t in ('phase_plus','phase_minus'))
    longchange={k:abs(long[k]/central[k]-1) for k in ('late_NX','late_EX','late_epsilon','late_NYL')}
    groups={
      'Jordan_geometry_and_gauge_algebra':dict(passed=max(errors)<1e-10 and max(wards)<1e-12,random_states=50,max_linear_system_error=max(errors),max_gauge_Ward_error=max(wards)),
      'reheating_conservation_and_radiation_reference':dict(passed=checks['max_Friedmann']<1e-6 and checks['max_Gauss']<1e-9 and checks['max_total_work']<1e-8 and checks['max_hidden_work']<1e-8 and max(rd_errors.values())<2e-4,checks=checks,RD_reference_relative_errors=rd_errors),
      'independent_Radau_tolerance_and_phase_controls':dict(passed=independent<2e-7 and tolerr<2e-7 and phase_response<1e-5,Radau_state_error=independent,tighter_relative_change=tolerr,phase_response_for_1e8_perturbation=phase_response,relative_charge_derivative_wrt_phase=phase_slope),
      'long_window_and_matter_growth_controls':dict(passed=longchange['late_NX']<.005 and longchange['late_epsilon']<.005 and growth_error<1e-8,long_relative_changes=longchange,matter_origin_growth_relative_error=growth_error,
          scope='Instantaneous vacuum criteria and homogeneous response only; finite-k/gauge/metric fluctuation certification remains open.')
    }
    write('REHEATING_VERIFICATION.json',dict(groups=groups,physical_modes_stability_certified=False))
    # Compact review data: dense early segment plus sparse late diagnostics.
    selected={}
    for tag in ('b2_g0.1','b10_g0.1','xi_zero_control','RD_reference'):
        arr=np.load(WORK/(tag+'.npz'))['samples'];indices=np.unique(np.r_[np.flatnonzero(arr[:,0]<=10)[::2],np.linspace(0,len(arr)-1,250).astype(int)])
        selected[tag]=arr[indices,:len(COLS)]
    np.savez_compressed(ROOT/'REHEATING_PROFILES.npz',**selected)
    print('REHEATING_VERIFY',(ROOT/'REHEATING_VERIFICATION.json').read_text(),flush=True)
    if not all(g['passed'] for g in groups.values()):raise SystemExit(1)


if __name__=='__main__':main()
