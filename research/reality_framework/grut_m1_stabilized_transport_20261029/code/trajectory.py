"""Build and run the scalar background and canonical transfer calculation."""
from pathlib import Path
import argparse,json,subprocess,time,hashlib
import numpy as np
from scipy.integrate._ivp import dop853_coefficients as dc
import model as m
ROOT=m.ROOT;W=ROOT/'.work';W.mkdir(exist_ok=True)

def build():
 vals={'LX':m.lx,'LY':m.ly,'LS':m.ls,'MM':m.M,'HH':m.h,'MU':m.mu,'KAP':m.k,'XI':m.xi}
 (W/'values.h').write_text(''.join(f'#define {k} ({v:.17g})\n' for k,v in vals.items()))
 def fmt(n,a):
  a=np.asarray(a);shape=''.join(f'[{x}]' for x in a.shape)
  v=','.join(f'{x:.17g}' for x in a) if a.ndim==1 else ','.join('{'+','.join(f'{x:.17g}' for x in row)+'}' for row in a)
  return f'static const double {n}{shape}={{'+v+'};\n'
 (W/'dop853.h').write_text(''.join(fmt(n,a) for n,a in [('AA',dc.A[:12,:12]),('BB',dc.B),('CC',dc.C[:12]),('EE3',dc.E3),('EE5',dc.E5)]))
 p=subprocess.run(['gcc','-O3','-std=c11','-Wall','-Wextra','-I',str(W),str(ROOT/'code/evolve.c'),'-lm','-o',str(W/'evolve')],capture_output=True,text=True)
 (W/'compile.log').write_text(p.stdout+p.stderr)
 if p.returncode:raise RuntimeError(p.stderr)
 return W/'evolve'

def run(tag,end=800.,rtol=2e-9,ks=(0.,),eg=m.e,lock=m.LOCK,z=None,delta=None):
 exe=build();v=np.zeros(17);v[0]=np.sqrt(6)/3;v[4:10]=m.real(m.initial() if z is None else z)
 if delta is not None:v[4:16]+=delta
 times=np.unique(np.r_[np.linspace(0,min(10,end),201),np.linspace(min(10,end),end,1601)])
 with (W/(tag+'.in')).open('wb') as f:
  np.array([len(ks),len(times)],dtype='<i4').tofile(f);np.array([end,rtol,eg,lock.real,lock.imag,.08/3],dtype='<f8').tofile(f);v.astype('<f8').tofile(f);np.array(ks,dtype='<f8').tofile(f);times.astype('<f8').tofile(f)
 tic=time.monotonic();p=subprocess.run([str(exe),str(W/(tag+'.in')),str(W/(tag+'.bin'))],capture_output=True,text=True)
 (W/(tag+'.log')).write_text(p.stdout+p.stderr)
 if p.returncode:raise RuntimeError(f'{tag} status {p.returncode} '+p.stderr)
 a=np.fromfile(W/(tag+'.bin'),dtype='<f8').reshape(len(times),18+144*len(ks));bg=a[:,1:18];trans=a[:,18:].reshape(len(times),len(ks),12,12)
 obs=[]
 for t,v in zip(times,bg):
  I,J,rad,lna=v[:4];aa=np.exp(lna);H=np.sqrt((J*J/2+4.5*I*I+rad)/3);R=18*I*I-J*J;z=m.complexvec(v[4:10]);vel=m.complexvec(v[10:16]);V,F,VR=m.evaluate(z,R,lock);charges=-2*aa**3*np.imag(np.conj(z)*vel);E=.5*(v[10:16]@v[10:16])+V
  obs.append([t,aa,H,R,charges[0],charges@[1,3,6],E,E-v[16],1-2*m.eta**2*VR])
 obs=np.array(obs);late=times>=.8*end;late2=times>=.9*end;mean=lambda mask:float(np.trapezoid(obs[mask,4],times[mask])/(times[mask][-1]-times[mask][0]))
 out={'tag':tag,'end_tHi':end,'rtol':rtol,'e':eg,'ks_over_Hi':list(ks),'lock':[lock.real,lock.imag],
 'late_mean_NX':mean(late),'late_half_window_mean_NX':mean(late2),'max_abs_gauss':float(max(abs(obs[:,5]))),
 'relative_energy_ledger_error':float(max(abs(obs[:,7]-obs[0,7]))/max(abs(obs[0,6]),1)),
 'min_effective_Planck_ratio':float(min(obs[:,8])),'max_effective_Planck_ratio':float(max(obs[:,8])),
 'Ricci_range': [float(min(obs[:,3])),float(max(obs[:,3]))],'seconds':time.monotonic()-tic,'engine':p.stdout.strip(),
 'scope':'Prescribed reheating background and scalar/longitudinal tangent. No metric, inflaton or radiation perturbations; no full CDI transfer or state law.'}
 np.savez_compressed(ROOT/(tag+'.npz'),times=times,bg=bg,trans=trans,obs=obs,ks=np.array(ks));(ROOT/(tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2),flush=True);return out

def main():
 p=argparse.ArgumentParser();p.add_argument('--tag',default='holder');p.add_argument('--end',type=float,default=800);p.add_argument('--rtol',type=float,default=2e-9);p.add_argument('--ks',default='0');p.add_argument('--original',action='store_true');a=p.parse_args()
 z=np.array(m.P['preparation']['initial_fields_H_units']);z=m.complexvec(z*np.sqrt(2)) if a.original else None
 ks=[] if a.ks=='none' else [float(x) for x in a.ks.split(',')]
 run(a.tag,a.end,a.rtol,ks,lock=0j if a.original else m.LOCK,z=z)
if __name__=='__main__':main()
