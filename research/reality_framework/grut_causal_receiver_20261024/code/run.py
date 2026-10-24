#!/usr/bin/env python3
"""Execute causal-layer free packets and unchanged prepared radial collapse."""
import argparse,hashlib,json,subprocess,tempfile,time
from pathlib import Path
import numpy as np
from scipy.fft import dst
from scipy.integrate import solve_bvp

ROOT=Path(__file__).resolve().parents[1]
BASE=dict(h=.25,dt=.03125,length=450.,tmax=12000.,probe=240.,absorber=360.,sigma=20.,alpha=.1,nonlinear=True,mass=1.)
CASES={
 'P025':{},'P0125':dict(h=.125,dt=.015625),'P025_time':dict(dt=.015625),
 'P00625':dict(h=.0625,dt=.015625),
 'P025_remote':dict(length=810.,absorber=720.,tmax=24000.),
 'P025_remote_time':dict(length=810.,absorber=720.,tmax=12000.,dt=.015625),
 'free_X':dict(h=.5,dt=.03125,tmax=6000.,nonlinear=False,k=.1,packet_width=40.,center=120.,reference_length=1800.),
 'free_X_time':dict(h=.5,dt=.015625,tmax=6000.,nonlinear=False,k=.1,packet_width=40.,center=120.,reference_length=1800.),
 'free_s':dict(h=.5,dt=.015625,tmax=6500.,nonlinear=False,mass=np.sqrt(6.),k=.2,packet_width=30.,center=120.,reference_length=1800.),
 'free_soft':dict(h=2.,dt=.0625,length=1090.,absorber=1000.,probe=800.,tmax=60000.,nonlinear=False,k=.02,packet_width=100.,center=300.,reference_length=5000.),
 'smoke':dict(h=.5,dt=.03125,tmax=100.),
}


def transform(z):return dst(z.real,type=1,norm='ortho')+1j*dst(z.imag,type=1,norm='ortho')


def dispersion(k,h,mass):
    sine=np.sin(k*h/2);p=2*sine/h*(1+sine*sine/6)
    return np.sqrt(mass*mass+p*p)


def free_initial(c):
    h=c['h'];r=np.arange(1,round(c['reference_length']/h))*h
    center,width,k=c['center'],c['packet_width'],c['k']
    u=np.exp(-(r-center)**2/(2*width*width)+1j*k*(r-center))-np.exp(-(r+center)**2/(2*width*width)-1j*k*(r+center))
    modes=transform(u);ks=np.pi*np.arange(1,len(r)+1)/c['reference_length'];omega=dispersion(ks,h,c['mass'])
    v=transform(-1j*omega*modes)
    count=round(c['length']/h)-1
    return r[:count],u[:count],v[:count],dict(r=r,modes=modes,omega=omega,k=ks)


def initial(c,path):
    h,L=c['h'],c['length'];r=np.arange(1,round(L/h))*h
    if c['nonlinear']:
        f=lambda x:.18*np.exp(-x*x/(2*60.**2))
        rhs=lambda x,y:np.array([y[1],(3*(y[0]**2-1)+f(x)**2)*y[0]])
        x=np.linspace(0,L,601);guess=np.array([np.sqrt(1-f(x)**2/3),f(x)**2*x/(3*60.**2)])
        solution=solve_bvp(rhs,lambda left,right:[left[1],right[0]-1],x,guess,S=np.diag([0.,-2.]),tol=1e-10,bc_tol=1e-11,max_nodes=12000)
        assert solution.success
        field=f(r);u=r*field;v=-1j*u;neutral=solution.sol(r)[0];reference=None
    else:
        r,u,v,reference=free_initial(c);neutral=np.ones_like(r)
    table=np.column_stack([r,u.real/r,u.imag/r,neutral,v.real/r,v.imag/r,np.zeros_like(r)])
    np.savetxt(path,table,delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
    charge=4*np.pi*h*np.sum(-np.imag(u.conj()*v))
    outside=r>=c['absorber'];tail=4*np.pi*h*np.sum(-np.imag(u[outside].conj()*v[outside]))
    return dict(initial_numeric_sha256=hashlib.sha256(table.tobytes()).hexdigest(),
        Qbar_initial_from_data=float(charge),initial_charge_fraction_in_layer=float(tail/charge),
        initial_memory='psi=q=0 in the direct second-order causal formulation. Initial nonzero physical tails are retained and measured, not removed.'),reference


def gradient(z,h):
    u=np.r_[0j,z,0j];extended=np.r_[-u[1],u,-u[-2]]
    return ((9/8)*(u[1:]-u[:-1])-(1/24)*(extended[3:]-extended[:-3]))/h


def free_comparison(c,reference,final,initial_E):
    time_final=c['tmax'];dt=c['dt'];steps=round(time_final/dt)
    omega=reference['omega'];modes=reference['modes'];n=len(final);h=c['h']
    actual_u=final[:,0]*(final[:,1]+1j*final[:,2]);actual_v=final[:,0]*(final[:,4]+1j*final[:,5])
    x=-1j*omega*dt;z=1+x+x*x/2+x**3/6+x**4/24
    mask_node=final[:,0]<c['absorber']-4*h
    mask_face=(np.arange(n+1)+.5)*h<c['absorber']-4*h
    comparisons={}
    for name,evolution in [('exact_time',np.exp(-1j*omega*time_final)),('same_RK4_bulk',z**steps)]:
        exact_u=transform(modes*evolution)[:n];exact_v=transform(-1j*omega*modes*evolution)[:n]
        du=actual_u-exact_u;dv=actual_v-exact_v
        error=2*np.pi*h*(sum((abs(dv)**2+c['mass']**2*abs(du)**2)[mask_node])+sum((abs(gradient(du,h))**2)[mask_face]))
        comparisons[name+'_physical_energy_error_over_initial']=float(error/initial_E)
    ks=reference['k'];sine=np.sin(ks*h/2);p=2*sine/h*(1+sine*sine/6)
    velocity=p/omega*np.cos(ks*h/2)*(1+.5*sine*sine)
    fast=velocity*time_final+c['center']+6*c['packet_width']>c['reference_length']
    spectral_energy=omega**2*abs(modes)**2
    comparisons['reference_fast_spectral_tail_energy_fraction']=float(sum(spectral_energy[fast])/sum(spectral_energy))
    comparisons['scope']='Physical-region error against a large-domain sine-mode solution. Same-RK4 comparison separates the boundary from the bulk time error. Very fast reference-tail contamination is estimated separately.'
    return comparisons


def run(name,binary):
    c={**BASE,**CASES[name]};data=ROOT/'data';data.mkdir(exist_ok=True);prefix=data/name
    with tempfile.TemporaryDirectory(prefix='grut-cfs-initial-') as tmp:
        init=Path(tmp)/'initial.csv';preparation,reference=initial(c,init)
        command=[str(binary),str(init),str(prefix)]+[str(c[k]) for k in ['h','dt','length','tmax','probe','absorber','sigma','alpha']]+[str(int(c['nonlinear'])),str(c['mass'])]
        print('START',name,c,flush=True);start=time.monotonic()
        with (data/f'{name}.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    history=np.loadtxt(str(prefix)+'_history.csv',delimiter=',',skiprows=1)
    trace=np.loadtxt(str(prefix)+'_trace.csv',delimiter=',',skiprows=1)
    final=np.loadtxt(str(prefix)+'_final.csv',delimiter=',',skiprows=1)
    memory=np.loadtxt(str(prefix)+'_memory.csv',delimiter=',',skiprows=1)
    hashes={}
    for key,array in [('trace',trace),('final',final),('memory',memory),('history',history)]:
        hashes[key]=hashlib.sha256(array.tobytes()).hexdigest()
        np.savez_compressed(str(prefix)+f'_{key}.npz',**{key:array})
    for key in ['trace','final','memory','history']:Path(str(prefix)+f'_{key}.csv').unlink()
    e0,q0=history[0,1:3];late=history[:,0]>.8*c['tmax'];keep=np.unique(np.r_[np.arange(0,len(history),4),len(history)-1])
    out=dict(name=name,parameters=c,seconds=time.monotonic()-start,preparation=preparation,
        source_sha256=hashlib.sha256((ROOT/'code/radial_cfs.cpp').read_bytes()).hexdigest(),
        raw_numeric_sha256=hashes,Ebar_initial=float(e0),Qbar_initial=float(q0),
        maximum_energy_plus_receiver_relative_error=float(max(abs((history[:,1]+history[:,3])/e0-1))),
        maximum_charge_plus_receiver_relative_error=float(max(abs((history[:,2]+history[:,4])/q0-1))),
        core60_late_mean_fraction=float(np.mean(history[late,5])/q0),core60_final_fraction=float(history[-1,5]/q0),
        maximum_layer_charged_mass2_shift=float(max(history[:,6])),
        history_columns=['tau','E','Q','E_receiver','Q_receiver','Q60','max_layer_mass2_shift','Fmax'],
        history=history[keep].tolist(),
        scope='Prepared spherical CFS evolution, with signed numerical boundary work. The receiver ledger alone is not a positive particle spectrum or actual cosmological source.')
    if reference is not None:out['free_comparison']=free_comparison(c,reference,final,e0)
    (ROOT/f'{name}_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('DONE',name,'seconds',out['seconds'],'Eerr',out['maximum_energy_plus_receiver_relative_error'],
        'Qerr',out['maximum_charge_plus_receiver_relative_error'],'retention',out['core60_late_mean_fraction'],out.get('free_comparison',{}),flush=True)
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('names',nargs='+',choices=list(CASES));args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='grut-cfs-build-') as tmp:
        binary=Path(tmp)/'radial';subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'code/radial_cfs.cpp'),'-o',str(binary)],check=True)
        for name in args.names:run(name,binary)


if __name__=='__main__':main()
