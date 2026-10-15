"""Full-law squared-action integrals with quadrature and state-error intervals.

The midpoint expansion integrates f=||L U_S(t)v||^2 through fourth order.
A sixth-derivative bound controls its remainder. Independent outgoing norm
upper sums bound the full-H state error and its time-integrated square.
"""
from math import factorial,ceil,comb,sqrt
import json,time
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply
from boundary import RF
import sys
sys.path.insert(0,str(RF/'cp12_shortcut_dynamics_20261015'))
from integrals import schur_norm

ALLOWANCE=5e-10


def dot(a,b):return np.vdot(a,b).real


class ActionEvolution:
    def __init__(self,H,B,operators,diagonals,diagonal_bounds,q,shift,dt=1/32,batch=2,hopping_norm_bound=None,matrix_free=False):
        if matrix_free:
            from linear_propagator import shifted_operators
            self.A,self.G,self.trace=shifted_operators(H,shift)
        else:
            self.A=H-shift*eye(H.shape[0],format='csr');self.G=-1j*self.A
            self.trace=-1j*complex(self.A.diagonal().sum())
        self.B=B;self.bn=schur_norm(B)
        self.operators=operators;self.diagonals=diagonals;self.dt=dt;self.batch=batch
        self.hopping_norm_bound=2*q if hopping_norm_bound is None else hopping_norm_bound
        self.bounds={name:self.hopping_norm_bound**2 for name in operators}|diagonal_bounds

    def run(self,v,times,initial_error=0.):
        assert v.ndim==1
        moment=[float(np.linalg.norm(v))];w=v.copy()
        for k in range(1,7):w=self.A@w;moment.append(float(np.linalg.norm(w)))
        sixth=sum(comb(6,k)*moment[k]*moment[6-k] for k in range(7))
        integrals={key:0. for key in self.bounds};quad={key:0. for key in self.bounds}
        epsilon5=epsilon1=error_l2=0.;clock=0.;out={};last=time.monotonic()
        for end in times:
            steps=max(1,ceil((end-clock)/self.dt));h=(end-clock)/steps
            coefficients=[h**(k+1)/(2**k*factorial(k+1)) for k in range(6)]
            tail=h**7/(2**6*factorial(7))*self.bn*moment[6]
            for j in range(0,steps,self.batch):
                count=min(self.batch,steps-j)
                grid=expm_multiply(self.G,v,start=0,stop=count*h,num=2*count+1,traceA=self.trace)
                for k in range(count):
                    ws=[grid[2*k+1]]
                    for order in range(1,6):ws.append(self.A@ws[-1])
                    local=low=0.
                    for order,c in enumerate(coefficients):
                        x=c*float(np.linalg.norm(self.B@ws[order]));local+=x
                        if order<2:low+=x
                    epsilon5+=local+tail;epsilon1+=low+h**3/24*self.bn*moment[2]
                    error=initial_error+min(epsilon5,epsilon1)
                    error_l2+=h*error*error
                    for name,blocks in self.operators.items():
                        f0=f2=f4=0.
                        for L in blocks:
                            z=[L@ws[a] for a in range(5)]
                            f0+=dot(z[0],z[0]);f2+=2*(dot(z[1],z[1])-dot(z[0],z[2]))
                            f4+=2*dot(z[0],z[4])-8*dot(z[1],z[3])+6*dot(z[2],z[2])
                        integrals[name]+=h*f0+h**3/24*f2+h**5/1920*f4
                    for name,D in self.diagonals.items():
                        f0=dot(ws[0],D*ws[0]);f2=2*(dot(ws[1],D*ws[1])-dot(ws[0],D*ws[2]))
                        f4=2*dot(ws[0],D*ws[4])-8*dot(ws[1],D*ws[3])+6*dot(ws[2],D*ws[2])
                        integrals[name]+=h*f0+h**3/24*f2+h**5/1920*f4
                    for name in quad:quad[name]+=h**7/(2**6*factorial(7))*self.bounds[name]*sixth
                v=grid[-1].copy()
                if time.monotonic()-last>15:
                    print(json.dumps(dict(stage='action_integrals',time=clock+(j+count)*h,
                        dimension=len(v),state_error=error,outside=integrals['outside'],inside=integrals['inside'])),flush=True);last=time.monotonic()
            clock=end;values={}
            for name,central in integrals.items():
                qerr=quad[name]+ALLOWANCE
                retained=[max(0.,central-qerr),max(0.,central+qerr)]
                radius=sqrt(self.bounds[name]*error_l2)+ALLOWANCE
                values[name]=dict(value=central,quadrature_error=quad[name],
                    retained_quadrature_interval=retained,
                    full_law_interval=[max(0.,sqrt(retained[0])-radius)**2,(sqrt(retained[1])+radius)**2],
                    action_L2_error_radius=radius)
            error=initial_error+min(epsilon5,epsilon1)
            final={name:dict(value=float(dot(v,D*v)),
                interval=[max(0.,max(0.,float(dot(v,D*v)))**.5-sqrt(self.bounds[name])*error-ALLOWANCE)**2,
                    (max(0.,float(dot(v,D*v)))**.5+sqrt(self.bounds[name])*error+ALLOWANCE)**2])
                for name,D in self.diagonals.items()}
            out[end]=dict(integrals=values,state=v.copy(),state_error=error,
                lag_state_error=min(epsilon5,epsilon1),integrated_state_error_square_upper=error_l2,
                final=final)
        return out
