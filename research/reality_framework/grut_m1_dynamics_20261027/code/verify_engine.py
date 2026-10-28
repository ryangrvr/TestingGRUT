"""Checks of the changed Hamiltonian/force and actual compiled evolution."""
from common import *
from scipy.integrate import solve_ivp

def evolve(tag,q,p,L,dt,T,sample=.04):
    d=WORK/tag;d.mkdir(exist_ok=True);initial=d/'initial.bin';prefix=d/'state'
    write_initial(initial,q,p,L)
    subprocess.run([str(build()),str(initial),str(prefix),str(dt),str(T),str(sample),str(T),'4','1','64'],check=True,stdout=subprocess.DEVNULL)
    data=np.genfromtxt(str(prefix)+'.tsv',names=True,delimiter='\t')
    files=sorted(d.glob('state_s*.bin'));N,_,_,a=snapshot(files[-1])
    return data,np.asarray(a).copy()

def main():
    rng=np.random.default_rng(40571);q=rng.normal(size=(3,7,7,7));q[2]+=7;p=rng.normal(size=q.shape)
    direction=rng.normal(size=q.shape);h=1e-5
    numerical=np.sum(potential(q+h*direction)-potential(q-h*direction))/(2*h)
    analytic=np.sum(potential_gradient(q)*direction)
    derivative_error=abs(numerical/analytic-1)
    L=13.;dx=L/7
    energy=totals(q,p,L)['E']
    integration_by_parts=float((.5*np.sum(p*p)+potential(q).sum()-.5*np.sum(q*lap(q,L)))*dx**3)
    spatial_error=abs(energy/integration_by_parts-1)
    # Independent 0D DOP853 orbit, including the new cubic neutral force.
    q0=np.array([13.,.2,6.]);p0=np.array([.3,9.,-.1]);N=6
    qa=np.broadcast_to(q0[:,None,None,None],(3,N,N,N)).copy()
    pa=np.broadcast_to(p0[:,None,None,None],qa.shape).copy()
    def rhs(t,y):return np.r_[y[3:],-potential_gradient(y[:3])]
    exact=solve_ivp(rhs,(0,2),np.r_[q0,p0],method='DOP853',rtol=2e-13,atol=2e-13).y[:,-1]
    errors=[];cpp_energy_errors=[];drifts=[]
    for dt in (.04,.02,.01):
        data,a=evolve('verify_uniform_'+str(dt),qa,pa,12.,dt,2.)
        got=a[:,0,0,0];errors.append(float(np.max(abs(got-exact))))
        cpp_energy_errors.append(float(abs(data['rho'][0]*12**3/totals(qa,pa,12.)['E']-1)))
        drifts.append(float(np.max(abs(data['nQ']/data['nQ'][0]-1))))
    order_ratios=[errors[i]/errors[i+1] for i in range(2)]
    # Independent growing eigenmode on a uniform rotating cold fluid.
    # Use the actual lattice k symbol, not a continuum frequency at finite spacing.
    N=24;L=48.;F=1.;S=brentq(lambda s:3*(s*s-1)*s-.6*(s-1)**2+F*F*s,0,1)
    om=np.sqrt(S*S+.1*F*F);k=2*np.pi*3/L;dx=L/N
    kd=4*np.sin(k*dx/2)**2*(1+np.sin(k*dx/2)**2/3)/dx**2
    D=3*(3*S*S-1)-1.2*(S-1)+F*F
    # y=(u,v,s,udot,vdot,sdot), all canonical amplitudes in the rotating frame.
    A=np.zeros((6,6));A[:3,3:]=np.eye(3)
    A[3,:]=[-kd-.2*F*F,0,-2*F*S,0,2*om,0]
    A[4,:]=[0,-kd,0,-2*om,0,0]
    A[5,:]=[-2*F*S,0,-kd-D,0,0,0]
    ev,vec=np.linalg.eig(A);idx=np.argmax(ev.real);gamma=float(ev[idx].real)
    z=vec[:,idx].real;z/=np.linalg.norm(z[:3]);amp=1e-5
    mode=np.cos(k*np.arange(N)*dx)[:,None,None];mode=np.broadcast_to(mode,(N,N,N))
    q=np.array([np.full_like(mode,V*F),np.zeros_like(mode),np.full_like(mode,V*S)])+amp*z[:3,None,None,None]*mode
    p=np.array([np.zeros_like(mode),np.full_like(mode,om*V*F),np.zeros_like(mode)])
    p+=amp*np.array([z[3]-om*z[1],z[4]+om*z[0],z[5]])[:,None,None,None]*mode
    T=16.;data,a=evolve('verify_growing_mode',q,p,L,.02,T,sample=.1)
    qt=a[:3];pt=a[3:];c=np.cos(om*T);s=np.sin(om*T)
    rotating=np.array([c*qt[0]+s*qt[1]-V*F,-s*qt[0]+c*qt[1],qt[2]-V*S])
    projection=2*np.mean(rotating*mode,axis=(1,2,3))
    predicted=amp*z[:3]*np.exp(gamma*T)
    growth_error=float(np.linalg.norm(projection-predicted)/np.linalg.norm(predicted))
    out={'potential_directional_derivative_relative_error':derivative_error,
         'spatial_Hamiltonian_integration_by_parts_error':spatial_error,
         'independent_DOP853_uniform_orbit_errors':errors,'time_order_ratios':order_ratios,
         'compiled_vs_python_initial_energy_errors':cpp_energy_errors,
         'uniform_charge_drifts':drifts,
         'lattice_cold_mode':{'k':k,'k_lattice_squared':kd,'predicted_growth':gamma,
                              'evolved_time':T,'projection_relative_error':growth_error},
         'pass':bool(derivative_error<1e-7 and spatial_error<1e-12 and max(cpp_energy_errors)<1e-12
                     and 14<min(order_ratios)<18 and errors[-1]<1e-6 and growth_error<.002)}
    save('ENGINE_CHECKS.json',out);print(json.dumps(out,indent=2))
    if not out['pass']:raise RuntimeError('Engine verification failed')

if __name__=='__main__':main()
