"""Run the declared Jordan-frame reheating models with the full tree sector."""
from pathlib import Path
import json,math,subprocess,hashlib
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from scipy.integrate._ivp import dop853_coefficients as dc
from common import ROOT,WORK,write,previous

INP=previous();P=INP['normalized_parameters']
LX,LY,LS=.001,1.,1.;M,HH,MU,KAP=30.,40.,INP['dimensionless_early']['mu'],P['kappa']
H0=P['H_GeV'];MRED=2.435e18;ETA2=(H0/MRED)**2;XI=-5/6
COLS=['tau','a','H','R','w_inflaton','radiation_fraction','NX','Nlight','Gauss','EX','epsilon','NYH','NYL','NSigma','Sigma_over_vev','Friedmann_residual','rho_minus_work','hidden_energy_minus_work','F_over_M2','trace_residual','curvature_mass2','fixed_Y_X_Hessian_min','rho_rad','rho_total','kappa_interaction','mu_interaction','absX','absY','hidden_canonical_energy']


def potential_force(z):
    x,y,s=z;xx,yy,ss=abs(z)**2
    V=LX*xx*xx+LY*yy*yy+M*M*yy+LS*(ss-HH*HH/2)**2+2*np.real(KAP*y.conjugate()*x**3+MU*s.conjugate()*y*y)
    F=np.array([2*LX*xx*x+3*KAP*y*x.conjugate()**2,
                (M*M+2*LY*yy)*y+KAP*x**3+2*MU*s*y.conjugate(),
                2*LS*(ss-HH*HH/2)*s+MU*y*y])
    return V,F


def geometric(v,beta,xi,eta2):
    z=v[2:8:2]+1j*v[3:8:2];vel=v[8:14:2]+1j*v[9:14:2]
    pot,force=potential_force(z);K=np.vdot(vel,vel).real
    f=1-2*xi*eta2*abs(z[0])**2;fd=-4*xi*eta2*np.real(z[0].conjugate()*vel[0])
    R=(2*beta**2*v[0]**2-v[1]**2+eta2*(4*pot-2*K+12*xi*abs(vel[0])**2-12*xi*np.real(z[0].conjugate()*force[0])))/(f+12*xi**2*eta2*abs(z[0])**2)
    ax=-3*v[16]*vel[0]-force[0]-xi*R*z[0]
    fdd=-4*xi*eta2*(abs(vel[0])**2+np.real(z[0].conjugate()*ax))
    hd=-(v[1]**2+2*eta2*K+4*v[14]/3+fdd-v[16]*fd)/(2*f)
    return z,vel,pot,force,f,fd,R,hd


def rhs(t,v,beta,decay,xi,eta2):
    z,vel,pot,F,f,fd,R,hd=geometric(v,beta,xi,eta2);out=np.zeros(19)
    out[0]=v[1];out[1]=-(3*v[16]+decay)*v[1]-beta**2*v[0]
    acc=-3*v[16]*vel-F;acc[0]-=xi*R*z[0]
    out[2:8:2]=vel.real;out[3:8:2]=vel.imag;out[8:14:2]=acc.real;out[9:14:2]=acc.imag
    K=np.vdot(vel,vel).real
    out[14]=decay*v[1]**2-4*v[16]*v[14];out[15]=v[16];out[16]=hd
    out[17]=-3*v[16]*(v[1]**2+2*eta2*K+4*v[14]/3)+.5*fd*R
    out[18]=-6*v[16]*K-2*xi*R*np.real(z[0].conjugate()*vel[0])
    return out


def prepared(beta,delta=0.,xi=XI,eta2=ETA2,rd=False):
    direction=-1j*np.exp(1j*delta)
    base=np.array([76.39978266076902,-2.111509854283515,-.9393290404896636,28.293142125868524])
    def state(q):
        z=np.array([q[0],q[1]+1j*q[2],direction*q[3]],complex);pot,_=potential_force(z)
        f=1-2*xi*eta2*abs(z[0])**2;v=np.zeros(19);v[2:8:2]=z.real;v[3:8:2]=z.imag;v[16]=1.
        if rd:v[14]=3*f-eta2*pot
        else:v[0]=math.sqrt(2*(3*f-eta2*pot))/beta
        return v
    if not rd:
        def eq(q):
            v=state(q);z,vel,pot,F,f,fd,R,hd=geometric(v,beta,xi,eta2)
            F=F.copy();F[0]+=xi*R*z[0]
            return [F[0].real,F[1].real,F[1].imag,np.real(F[2]*direction.conjugate())]
        sol=root(eq,base,tol=1e-12);base=sol.x;res=float(np.linalg.norm(eq(base)))
        assert res<1e-7,res
    else:res=0.
    return state(base),res


def build():
    b=WORK/'build';b.mkdir(exist_ok=True)
    def fmt(name,a):
        a=np.asarray(a)
        if a.ndim==1:return 'static const double '+name+'['+str(len(a))+']={'+','.join(format(x,'.17g') for x in a)+'};\n'
        return 'static const double '+name+'['+str(a.shape[0])+']['+str(a.shape[1])+']={'+','.join('{'+','.join(format(x,'.17g') for x in row)+'}' for row in a)+'};\n'
    (b/'dop853.h').write_text(''.join(fmt(n,a) for n,a in [('AA',dc.A[:12,:12]),('CC',dc.C[:12]),('BB',dc.B),('EE3',dc.E3),('EE5',dc.E5)]))
    exe=b/'reheat';subprocess.run(['gcc','-O3','-std=c11','-Wall','-Wextra','-I',str(b),str(ROOT/'code/reheat.c'),'-lm','-o',str(exe)],check=True)
    return exe


def run(exe,tag,beta,decay,xi=XI,eta2=ETA2,delta=0.,end=None,tol=2e-11,samples=10000,rd=False,frozen_prep=False):
    prep_xi=XI if frozen_prep else xi
    v,res=prepared(beta,delta,prep_xi,eta2,rd)
    if frozen_prep:
        z=v[2:8:2]+1j*v[3:8:2];pot,_=potential_force(z)
        v[0]=math.sqrt(2*(3*(1-2*xi*eta2*abs(z[0])**2)-eta2*pot))/beta
    end=end or 25/decay
    args=[beta,decay,xi,eta2,MU,KAP,end,tol,tol*.001]+v[2:8].tolist()+[v[0],v[14],samples]
    p=subprocess.run([str(exe),*[format(x,'.17g') for x in args]],capture_output=True,text=True,check=True)
    data=[];check=None
    for line in p.stdout.splitlines():
        q=line.split()
        if q[0]=='S':data.append([float(x) for x in q[1:]])
        elif q[0]=='CHECK':check=[float(x) for x in q[1:]]
    arr=np.array(data);np.savez_compressed(WORK/(tag+'.npz'),samples=arr)
    assert check is not None
    late=arr[:,0]>=.8*end
    def average(col,mask=late):return float(np.trapezoid(arr[mask,col],arr[mask,0])/(arr[mask,0][-1]-arr[mask,0][0]))
    a=arr[-1,1];rad=arr[-1,22]
    rho_rad=rad*MRED*MRED*H0*H0
    T=(30*rho_rad/(math.pi**2*106.75))**.25
    sentropy=2*math.pi**2/45*106.75*T**3
    nx=average(6);yieldQ=nx*H0**3/(a**3*sentropy)
    out=dict(tag=tag,beta=beta,Gamma_over_H=decay,xi_X=xi,eta_squared=eta2,delta=delta,rd_control=rd,end_tau=end,rtol=tol,samples=samples,
             prepared_conditional_residual=res,initial_state=v.tolist(),final_state=arr[-1,len(COLS):].tolist(),
             final_a=a,radiation_fraction=arr[-1,5],final_H_over_Hi=arr[-1,2],Tvis_final_GeV=T,
             late_NX=nx,late_EX=average(9),late_epsilon=average(10),late_NYH=average(11),late_NYL=average(12),late_NSigma=average(13),
             yield_over_YQ=yieldQ/INP['benchmark']['YQ'],
             late_half_window_relative_charge_change=abs(average(6,arr[:,0]>=.9*end)/nx-1),
             checks=dict(steps=int(check[0]),rejected=int(check[1]),max_Friedmann=check[2],max_Gauss=check[3],max_trace=check[4],max_total_work=check[5],max_hidden_work=check[6],min_F=check[7]),
             raw_output_sha256=hashlib.sha256(p.stdout.encode()).hexdigest())
    if rd:
        mask=(arr[:,1]>=30)&(arr[:,1]<=40+1e-9)
        aa=arr[mask,1]
        out['reference_window_30_40']={name:float(np.trapezoid(arr[mask,COLS.index(name)],aa)/(aa[-1]-aa[0])) for name in ('NX','EX','epsilon','NYH','NYL','NSigma')}
    return out,arr


def main():
    exe=build();jobs=[]
    for beta in (2.,10.):
        for dec in (.03,.1,.3):jobs.append((f'b{beta:g}_g{dec:g}',beta,dec,{}))
    jobs += [('RD_reference',2.,.1,dict(rd=True,eta2=0.,end=799.5)),
             ('xi_zero_control',2.,.1,dict(xi=0.,frozen_prep=True)),
             ('central_tight',2.,.1,dict(tol=2e-12,samples=6000)),
             ('phase_plus',2.,.1,dict(delta=1e-8,samples=6000)),
             ('phase_minus',2.,.1,dict(delta=-1e-8,samples=6000))]
    rows=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(run,exe,tag,beta,dec,**kw):tag for tag,beta,dec,kw in jobs}
        for f in as_completed(futures):
            out,arr=f.result();rows.append(out)
            write('REHEATING_RESULTS.json',dict(rows=sorted(rows,key=lambda r:r['tag']),columns=COLS,scope=__doc__))
            print('REHEAT',out['tag'],'eps',out['late_epsilon'],'Y/YQ',out['yield_over_YQ'],'NY',out['late_NYH']+out['late_NYL'],'checks',out['checks'],flush=True)


if __name__=='__main__':main()
