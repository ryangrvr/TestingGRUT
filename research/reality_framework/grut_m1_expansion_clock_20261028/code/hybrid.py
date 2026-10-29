"""Quasistatic evolving-Q and radiation expansion kinetic diagnostics.

The stationary E(Q) spline determines the exchange energies. Optional
``--profiled`` uses seven computed stationary scattering profiles with linear
interpolation in Q; the default keeps Q=200000 angular strengths as a hybrid
sensitivity test. Frequency generators are linearly interpolated and off-grid
daughter deposition causes diffusion. Both effects are reported and tested.
Neither mode evolves the excited core, spatial inhomogeneity or formation.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq
import expansion
import transport

ROOT=Path(__file__).resolve().parents[1]

def profiled_operator(e,omega,ks,w):
    """Same transition/deposition algorithm as parent, with a supplied W(k,Q)."""
    interp=PchipInterpolator(ks,w,axis=0,extrapolate=False)
    n=len(e); L=np.zeros((2*n,2*n)); omitted=np.zeros(2*n)
    def get(k):
        if k<ks[0]: return interp(ks[0])*k/ks[0]
        if k>ks[-1]: raise ValueError('Core kernel support exceeded')
        return interp(k)
    for species in (0,1):
        for i,energy in enumerate(e):
            col=species*n+i;p=np.sqrt(energy*energy-1)
            ea=energy if species==1 else energy-2*omega
            if ea<=1:continue
            weights=get(np.sqrt(ea*ea-1));weights=(weights+weights.T)/2
            exchange=np.pi*weights[0,1]/(p*energy)
            eh=energy+(omega if species==1 else -omega)
            neutral=np.pi*weights[2,species]/(p*energy) if eh>np.sqrt(6.) else 0.
            eo=energy+(2*omega if species==1 else -2*omega)
            if e[0]-1e-11<=eo<=e[-1]+1e-11:
                L[(1-species)*n:(2-species)*n,col]+=exchange*expansion.deposit(e,[eo],[1.])
                L[col,col]-=exchange
            else:omitted[col]+=exchange
            if neutral:
                d,fraction=expansion.daughters(e,eh);omitted[col]+=neutral*(1-fraction)
                neutral*=fraction
                L[:n,col]+=neutral*d;L[n:,col]+=neutral*d;L[col,col]-=neutral
    invariant=np.r_[e-omega,e+omega]
    return L,omitted,float(np.max(abs(invariant@L)))


class Hybrid(expansion.Expanding):
    def __init__(self,e,mode,angles,epsilon,spacing,profiled=False):
        super().__init__(e,mode,angles,epsilon)
        d=json.loads((ROOT/'FAMILY.json').read_text())
        rows=d['nodes']
        qs=np.array([r['Q'] for r in rows])
        self.mass=CubicHermiteSpline(qs,[r['energy_mX'] for r in rows],[r['Omega'] for r in rows],extrapolate=False)
        self.qrange=(qs[0],qs[-1])
        self.profiled=profiled
        if profiled:
            self.qrange=(186000.,204000.)
            qgrid=np.arange(186000.,204001.,3000.)
            tables=[];ks=None
            for q in qgrid:
                z=np.load(ROOT/f'core_Q{int(q)}.npz')
                if ks is None:ks=z['k']
                elif not np.array_equal(ks,z['k']):raise ValueError('mismatched scattering k grids')
                tables.append(z['W'])
            tables=np.array(tables)
            omin=float(self.mass(self.qrange[1],1));omax=float(self.mass(self.qrange[0],1))
            self.os=np.unique(np.r_[omin,np.arange(omin+spacing,omax,spacing),omax,expansion.OMEGA])
        else:
            self.os=np.unique(np.r_[np.arange(.8800,.8870+spacing/2,spacing),expansion.OMEGA])
        Ls=[];omits=[];errs=[]
        original=transport.OMEGA
        try:
            for om in self.os:
                if profiled:
                    q=brentq(lambda q:float(self.mass(q,1))-om,*self.qrange,xtol=1e-7)
                    i=min(np.searchsorted(qgrid,q,side='right')-1,len(qgrid)-2)
                    frac=(q-qgrid[i])/(qgrid[i+1]-qgrid[i])
                    weights=(1-frac)*tables[i]+frac*tables[i+1]
                    L,o,err=profiled_operator(e,om,ks,weights)
                else:
                    transport.OMEGA=float(om)
                    L,o,_,err=transport.core_operator(e)
                Ls.append(L);omits.append(o);errs.append(err)
        finally:
            transport.OMEGA=original
        self.Ls=np.array(Ls);self.omits=np.array(omits)
        self.build_error=float(max(errs))
        self.nq=np.r_[-np.ones(self.n),np.ones(self.n)]

    def terms(self,y):
        n=self.n
        q=expansion.QREF+self.nq@y[:2*n]
        if not self.qrange[0] <= q <= self.qrange[1]:
            raise ValueError(f'Q={q} outside stationary family {self.qrange}')
        omega=float(self.mass(q,1))
        if not self.os[0] <= omega <= self.os[-1]:
            raise ValueError(f'omega={omega} outside generator grid')
        i=np.searchsorted(self.os,omega,side='right')-1
        i=min(i,len(self.os)-2)
        t=(omega-self.os[i])/(self.os[i+1]-self.os[i])
        L=self.Ls[i]*(1-t)+self.Ls[i+1]*t
        lost=self.omits[i]*(1-t)+self.omits[i+1]*t
        dL=(self.Ls[i+1]-self.Ls[i])/(self.os[i+1]-self.os[i])
        return q,omega,L,lost,dL

    def rhs(self,x,y):
        n=self.n; v=y[:2*n]
        factor=np.exp(-x)/self.epsilon
        q,om,L,omitted,dL=self.terms(v)
        c,lost=self.model.rates(v)
        c+=(L-self.model.L)@v
        lost+=float((omitted-self.model.omit_core)@v)
        red=np.r_[self.R@v[:n],self.R@v[n:]]
        work=self.work@(v[:n]+v[n:])
        ccore=L@v
        defect=self.e@(ccore[:n]+ccore[n:])+om*(self.nq@ccore)
        return np.r_[factor*c+red,factor*lost,work,factor*defect]

    def jac(self,x,y):
        n=self.n; v=y[:2*n]; factor=np.exp(-x)/self.epsilon
        q,om,L,omitted,dL=self.terms(v)
        J=self.model.rates(v,True)+(L-self.model.L)
        J+=np.outer((dL@v)*float(self.mass(q,2)),self.nq)
        z=np.zeros((2*n+3,2*n+3));z[:2*n,:2*n]=factor*J
        z[:n,:n]+=self.Rd;z[n:2*n,n:2*n]+=self.Rd
        z[-2,:n]=self.work;z[-2,n:2*n]=self.work
        return z

    def observables_hybrid(self,x,y):
        n=self.n;v=y[:2*n]
        q,om,L,_,_=self.terms(v)
        a=float(np.exp(x)); nplus=float(v[:n].sum());nminus=float(v[n:].sum())
        Eg=float(self.e@(v[:n]+v[n:]));estat=float(self.mass(q))
        return {'a_over_a0':a,'tau':float((1-1/a)/self.epsilon),
                'Nplus':nplus,'Nminus':nminus,'Qcore':float(q),'Omega':om,
                'gas_energy_mX':Eg,'gas_kinetic_mX':Eg-nplus-nminus,
                'work_redshift_mX':float(y[-2]),'omitted_events':float(y[-3]),
                'interpolated_generator_energy_defect_mX':float(y[-1]),
                'total_stationary_plus_gas_plus_work_mX':estat+Eg+float(y[-2]),
                'stationary_cold_excess_fraction':(estat+nplus+nminus)/expansion.EREF-1,
                'floor_particles':float(v[0]+v[n]),'minimum_bin':float(v.min())}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--shift-bins',type=int,default=32)
    p.add_argument('--angles',type=int,default=8)
    p.add_argument('--spacing',type=float,default=.00025)
    p.add_argument('--amax',type=float,default=2.)
    p.add_argument('--source',choices=['parent','h0025','h0025early','h005','h005early','h01'],default='parent')
    p.add_argument('--rtol',type=float,default=2e-7)
    p.add_argument('--profiled',action='store_true')
    p.add_argument('--resume',default=None,help='Existing compatible saved run to extend in ln(a).')
    p.add_argument('--tag',required=True)
    args=p.parse_args();t0=time.monotonic()
    e=expansion.grid(0,8.,'translated',args.shift_bins)
    source,audit=expansion.initial(e) if args.source=='parent' else expansion.radial_source(e,args.source)
    eps=expansion.clock()['H0_mX2_over_ncore0']
    model=Hybrid(e,'full',args.angles,eps,args.spacing,args.profiled)
    print('built',args.tag,'nodes',len(e),'frequency generators',len(model.os),'seconds',round(time.monotonic()-t0,2),flush=True)
    start_x=0.;y0=np.r_[source,0.,0.,0.]; prior_x=prior_y=None
    if args.resume:
        before=json.loads((ROOT/(args.resume+'.json')).read_text())
        old=before['configuration']
        for key in ('shift_bins','angles','source','profiled','spacing'):
            if old.get(key)!=getattr(args,key):raise ValueError('Incompatible restart setting: '+key)
        z=np.load(ROOT/(args.resume+'.npz'));prior_x=z['x'];prior_y=z['Y']
        if not np.array_equal(e,z['E']):raise ValueError('Incompatible restart grid')
        start_x=float(prior_x[-1]);y0=prior_y[:,-1]
        if start_x>=np.log(args.amax):raise ValueError('Restart must extend the scale factor')
    vals=np.unique(np.r_[start_x,-np.log(1-np.array([.1,1,10,100,1000,10000,100000])*eps),np.log([1.01,1.1,2.,5.,args.amax])])
    vals=vals[(vals>=start_x-1e-13)&(vals<=np.log(args.amax)+1e-13)]
    result=solve_ivp(model.rhs,(start_x,np.log(args.amax)),y0,method='BDF',
                     jac=model.jac,rtol=args.rtol,atol=1e-9,t_eval=vals)
    if not result.success:raise RuntimeError(result.message)
    tx,yy=result.t,result.y
    if prior_x is not None:
        tx=np.r_[prior_x,tx[1:]];yy=np.hstack((prior_y,yy[:,1:]))
    rows=[model.observables_hybrid(x,y) for x,y in zip(tx,yy.T)]
    ledger=np.array([r['total_stationary_plus_gas_plus_work_mX']-r['interpolated_generator_energy_defect_mX'] for r in rows])
    scope=('Evolving Omega(Q) and seven computed stationary scattering profiles, interpolated in Q and omega; frequency-interpolated generators on a fixed grid with known projection diffusion.' if args.profiled else
           'Hybrid evolving Omega(Q), but fixed reference-profile W(k); frequency-interpolated generators on a fixed grid, with known projection diffusion.')
    out={'configuration':vars(args),'scope':scope+' No excitation release or physical source beyond selected windows.',
         'solver':{'nfev':result.nfev,'njev':result.njev,'nlu':result.nlu,'seconds':time.monotonic()-t0},
         'checks':{'max_individual_generator_energy_error':model.build_error,
                   'max_energy_minus_generator_defect_drift_mX':float(np.max(abs(ledger-ledger[0]))),
                   'minimum_saved_population':float(result.y[:2*len(e)].min()),
                   'initial_source_audit':audit},'snapshots':rows}
    np.savez_compressed(ROOT/(args.tag+'.npz'),E=e,x=tx,Y=yy)
    (ROOT/(args.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('done',args.tag,'seconds',round(time.monotonic()-t0,2),'final',rows[-1],flush=True)

if __name__=='__main__':main()
