"""Closed, explicitly conditional moment controls, with active-set ODEs.

These are not the missing full momentum Boltzmann solution. Initial kinetic
equilibrium, fast s, vacuum masses and the formation recipes are supplied.
The hard-capture recycling limiter is a labelled transport approximation.
"""
import json,time,sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import ks_common as K
from thermo_state import Inverter,TH,initial_equilibrium
from stage2_moments import Case,Solver,capture_kernel
from bootstrap import ROOT,WORK

def save(name,data):
    def cv(v):
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError(type(v).__name__)
    (ROOT/name).write_text(json.dumps(data,indent=2,allow_nan=False,default=cv)+'\n')

class Control:
    def __init__(self,R,formation='gas',recycle=False,catalysis=False,energy='finite',rtol=2e-7,xend=25,tstop=.002,method='Radau',closure_stop=True):
        self.R=R;self.formation=formation;self.recycle=recycle;self.rtol=rtol;self.xend=xend;self.tstop=tstop
        self.method=method;self.closure_stop=closure_stop
        rule='a' if formation=='gas' else 'b'
        self.case=Case('control',R,rule=rule,e_cap_variant='omega' if energy=='linear' else 'finite',catalysis=catalysis)
        if rule=='a':self.case.cap_p=self.case.cap_m=1-1e-6
        self.helper=Solver(self.case);self.B=self.helper.B;self.c=self.case
        self.inv=Inverter(self.c.cap_p,self.c.e_cap,self.c.cap_m,self.c.e_cm)
        tab=K.RAT['channels_22']['XX_XX']['sigma_table'];xx=np.array(tab['sqrt_s_over_mX'])
        yy=np.array(tab['transport_nores']);ok=yy>0
        self.transport=PchipInterpolator(np.log(xx[ok]),np.log(yy[ok]),extrapolate=True)
        self.start=time.time();self.calls=0;self.active=[False,False]

    def state(self,x,w):
        nq=K.N_Q0*np.exp(-3*x)/K.MX**3
        t=np.exp(w[0]);up=self.c.cap_p if self.active[0] else min(w[1],1-1e-6)
        um=self.c.cap_m if self.active[1] else min(w[2],1-1e-6)
        nc=w[1]*nq if self.active[0] else 0.;nb=w[2]*nq if self.active[1] else 0.
        _,J,raw=self.inv._eval(w[0],up,um,nc,nb,*self.active,(0.,0.,0.))
        X,Xb,S=raw['X'],raw['Xb'],raw['S'];rho=X['rho']+Xb['rho']+S['rho']+self.c.e_cap*nc+self.c.e_cm*nb
        st=dict(t=t,up=up,um=um,us=up+um,n_c=nc,n_cb=nb,n_plus=X['n'],n_minus=Xb['n'],n_s=S['n'],
            rho_gas=X['rho']+Xb['rho']+S['rho'],rho_tot=rho,P_gas=X['P']+Xb['P']+S['P'],
            Q=X['n']+nc-Xb['n']-nb,P=Xb['n']+S['n']+nb,nq=nq)
        return st,J

    def recycling_rates(self,st):
        t=st['t'];sv,ein=capture_kernel(t,self.c.capture or 'hold');ep=ein+2*self.c.cap_p
        # Mean-energy transport estimate using the completed XX transport table.
        etarget=1+3*t;root_s=np.sqrt(2+2*ep*etarget)
        sigma=float(np.exp(self.transport(np.log(root_s))))
        gamma_soft=sigma*max(st['n_plus'],0)*K.MX**3
        pin=np.sqrt(max(ein*ein-1,0));pout=np.sqrt(ep*ep-1)
        reverse=sv*pin*ein/(pout*ep)
        gamma_rev=max(st['n_c'],0)*K.MX**3/K.Q_REF*reverse
        H=self.B.H_of_x(self.current_x)
        redshift=H*(ep*ep-1)/(ep*max(ein-1,1e-12))
        factor=(gamma_soft+redshift)/(gamma_soft+redshift+gamma_rev) if gamma_soft+redshift+gamma_rev>0 else 1.
        return factor,gamma_soft/H,gamma_rev/H

    def rhs(self,x,w):
        self.calls+=1;self.current_x=x;st,J=self.state(x,w);reactions=self.helper.rates(st)
        if self.calls%10000==0:
            print('PROGRESS',self.formation,self.recycle,self.calls,x,st['t'],self.active,st['n_c']/st['Q'],st['P']/st['Q'],flush=True)
        factor,el,rev=self.recycling_rates(st) if self.formation!='gas' else (1.,0.,0.)
        if self.recycle and 'capture' in reactions:
            g,dp,dnp=reactions['capture'];reactions['capture']=(factor*g,dp,dnp)
        H=self.B.H_of_x(x);source=sum(g*dp for g,dp,dn in reactions.values())/(H*K.MX**3)
        target=np.array([-3*st['Q'],-3*st['P']+source,-3*(st['rho_tot']+st['P_gas'])])
        for j,active in enumerate(self.active):
            if active:
                value=st['n_c'] if j==0 else st['n_cb']
                target+=3*value*J[:,j+1];J[:,j+1]*=st['nq']
        core=np.linalg.solve(J,target)
        return np.r_[core,3*np.exp(3*x)*st['P_gas']/self.E0]

    def jac(self,x,w):
        jac=np.zeros((4,4))
        for j in range(3):
            h=1e-5 if j==0 else (max(1e-12,abs(w[j])*1e-3) if self.active[j-1] else 1e-6)
            step=np.zeros(4);step[j]=h
            if j and self.active[j-1] and w[j]<=h:
                jac[:,j]=(self.rhs(x,w+step)-self.rhs(x,w))/h
            else:jac[:,j]=(self.rhs(x,w+step)-self.rhs(x,w-step))/(2*h)
        return jac

    def closure_rate(self,x,w):
        st,_=self.state(x,w);self.current_x=x
        if self.formation!='gas':return self.recycling_rates(st)[1]
        # Xbar relaxation on the more numerous X particles; a thermal
        # rate diagnostic, not a full nonzero-mu collision integral.
        f,lo,hi=self.helper.RATES.elastic['XXb_XXb'] if hasattr(self.helper,'RATES') else __import__('stage2_moments').RATES.elastic['XXb_XXb']
        sv=float(np.exp(f(np.clip(np.log(st['t']),np.log(lo),np.log(hi)))))
        return sv*st['n_plus']*K.MX**3/self.B.H_of_x(x)

    def run(self):
        rho=K.rho_onset(self.R)/K.MX**4;nq=K.N_Q0/K.MX**3
        initial=initial_equilibrium(nq,rho,self.c.cap_p,self.c.e_cap)
        t=initial['t'];up=initial['up'];um=initial['um'];nc=initial['n_cond']
        self.active=[nc>1e-12,False];w=np.array([np.log(t),nc/nq if self.active[0] else up,um,0.]);self.E0=rho
        x=0.;trajectories=[];switches=[];segments=[];stopped='xend'
        for segment in range(20):
            def cold(x,w):return w[0]-np.log(self.tstop)
            cold.direction=-1;cold.terminal=True
            def tiny(x,w):
                st,_=self.state(x,w);return st['P']/max(st['Q'],1e-300)-1e-18
            tiny.direction=-1;tiny.terminal=True
            events=[cold,tiny]
            for j in range(2):
                active=self.active[j];cap=self.c.cap_p if j==0 else self.c.cap_m
                def transition(x,w,j=j,active=active,cap=cap):return w[j+1] if active else w[j+1]-cap
                transition.terminal=True;transition.direction=-1 if active else 1;events.append(transition)
            def closure(x,w):return self.closure_rate(x,w)-1
            closure.direction=-1;closure.terminal=True
            if self.closure_stop:events.append(closure)
            start_active=list(self.active)
            sol=solve_ivp(self.rhs,(x,self.xend),w,method=self.method,jac=self.jac,rtol=self.rtol,atol=[1e-11,1e-15,1e-15,1e-11],
                          max_step=.08,events=events,dense_output=True)
            if not sol.success:raise RuntimeError(sol.message)
            samples=np.linspace(x,sol.t[-1],max(3,int((sol.t[-1]-x)*35)))
            for xx in samples:
                ww=sol.sol(xx);st,_=self.state(xx,ww);self.current_x=xx
                fac,el,rev=self.recycling_rates(st) if self.formation!='gas' else (1.,None,None)
                p=st['P']/st['Q'];energy=(np.exp(3*xx)*st['rho_tot']/self.E0+ww[3]-1)
                trajectories.append(dict(x=float(xx),T_over_m=st['t'],up=st['up'],um=st['um'],
                    pairs_per_net=p,f_at_stop=K.f_old(p),f_reference_Q200000=K.f_bench(p),
                    condensed_fraction=st['n_c']/st['Q'],anti_fraction=st['n_cb']/st['Q'],
                    charge_error=st['Q']/(nq*np.exp(-3*xx))-1,energy_work_error=energy,
                    recycling_factor=fac,hard_elastic_over_H=el,hard_reverse_over_H=rev,
                    energy_per_net_m=st['rho_tot']/st['Q'],kinetic_rate_over_H=self.closure_rate(xx,ww),
                    visible_T_GeV=self.B.T_of_x(xx)))
            segments.append(dict(start_x=x,end_x=float(sol.t[-1]),active=start_active,nfev=sol.nfev))
            x=float(sol.t[-1]);w=sol.y[:,-1]
            if sol.t_events[0].size:stopped=f'T={self.tstop}m';break
            if sol.t_events[1].size:stopped='pair_ratio=1e-18';break
            if self.closure_stop and sol.t_events[4].size:stopped='kinetic_shape_rate/H=1';break
            event=next((j for j in range(2) if sol.t_events[j+2].size),None)
            if event is None:break
            st,_=self.state(x,w)
            self.active[event]=not self.active[event]
            w[event+1]=0. if self.active[event] else (self.c.cap_p if event==0 else self.c.cap_m)
            switches.append(dict(x=x,species=event,new_active=self.active[event]))
        else:raise RuntimeError('Too many active-set transitions')
        final=trajectories[-1]
        emin=min(self.c.e_cap,self.c.e_cm,1.)
        final['conditional_energy_upper_f']=K.f_old(max(0.,(final['energy_per_net_m']-emin)/(2*emin)))
        result=dict(R=self.R,formation=self.formation,recycling_limiter=self.recycle,catalysis=self.case.catalysis,
            energy_variant=self.case.e_cap_variant,rtol=self.rtol,final=final,segments=segments,switches=switches,stop=stopped,
            checks=dict(max_charge_error=max(abs(r['charge_error']) for r in trajectories),
                        max_energy_work_error=max(abs(r['energy_work_error']) for r in trajectories)),
            trajectories=trajectories,calls=self.calls,seconds=time.time()-self.start,
            certified_full_kinetic_f_sym=None,
            scope='Conditional vacuum-mass fast-s moment model on the supplied radiation background; no full hard momentum network, physical formation kernel, in-medium EOS or certified kinetic threshold.')
        return result

if __name__=='__main__':
    formation=sys.argv[1] if len(sys.argv)>1 else 'instant';recycle='--recycle' in sys.argv
    R=float(sys.argv[2]) if len(sys.argv)>2 else K.R_VEC['e1e-7_band']+K.R_MED['e1e-7']
    name=f'CONTROL_DEEP_{formation}_{R:.5f}_recycle{int(recycle)}.json'
    try:
        out=Control(R,formation=formation,recycle=recycle).run();save(name,out);print(name,out['final'],out['checks'],out['seconds'],flush=True)
    except Exception as exc:
        save(name,dict(R=R,formation=formation,error=repr(exc),complete=False));raise
