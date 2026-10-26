"""Binding-field linear modes, conditional Gaussian seeds; no backreaction."""
import json,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.integrate import solve_ivp
from mode_common import ROOT,WORK,write
from modes import build,prepared,XI,ETA2,MU,KAP,geometric,matrices,rhs

def mass(v,b):
    _,_,_,_,_,_,R,_=geometric(v,b,XI,ETA2)
    return np.exp(2*v[15])*(.02*(v[2]**2+v[3]**2)-.03*(.924089637375362/133349619398.33566)**2-R/6)

def run(exe,r,ks,times,rtol=2e-10,suffix=''):
    tag=r['tag']+suffix;v,_=prepared(r['beta']);inp=WORK/f'phi_{tag}.in';out=WORK/f'phi_{tag}.bin'
    pars=[r['beta'],r['Gamma_over_H'],1e-7,XI,ETA2,MU,KAP,times[-1],rtol]
    with inp.open('wb') as f:
        np.array([len(ks),len(times)],dtype='<i4').tofile(f);np.array(pars,dtype='<f8').tofile(f)
        np.asarray(v,dtype='<f8').tofile(f);np.asarray(ks,dtype='<f8').tofile(f);np.asarray(times,dtype='<f8').tofile(f)
    proc=subprocess.run([str(exe),str(inp),str(out)],capture_output=True,text=True,check=True)
    raw=np.fromfile(out).reshape(len(times),20+4*len(ks));return raw[:,1:20],raw[:,20:].reshape(len(times),len(ks),2,2),proc.stdout.strip()

def main():
    base=json.loads((ROOT/'MODE_RESULTS.json').read_text());times=np.array(base['times'])
    ks=np.unique(np.r_[np.array(base['ks']),np.geomspace(.8,20,97)])
    exe=build();neutral=exe.parent/'neutral'
    subprocess.run(['gcc','-O3','-std=c11','-DNEUTRAL_SECTOR','-I',str(exe.parent),str(ROOT/'code/modes.c'),'-lm','-o',str(neutral)],check=True)
    rows=[];profiles=dict(times=times,ks=ks);maxsym=0.
    with ThreadPoolExecutor(max_workers=3) as pool:
        tasks={pool.submit(run,neutral,r,ks,times):r for r in base['rows']}
        for future in as_completed(tasks):
            r=tasks[future];bg,TT,check=future.result();m0=mass(bg[0],r['beta']);w0=np.sqrt(ks*ks+m0)
            # omega'_0=M_phi,conformal(0)/omega_0 since velocities/R'=0
            d0=m0/(2*w0*w0)
            covs=[np.array([[[1/(2*w),0],[0,w/2]] for w in w0]),
                  np.array([[[1/(2*w),-d/(2*w)],[-d/(2*w),(w*w+d*d)/(2*w)]] for w,d in zip(w0,d0)])]
            result=np.zeros((len(times),len(ks),3));ratios=[]
            for n,v in enumerate(bg):
                w2=ks*ks+mass(v,r['beta']);assert np.all(w2>0)
                w=np.sqrt(w2);EX=np.exp(4*v[15])*((v[8]**2+v[9]**2)+.001*(v[2]**2+v[3]**2)**2)
                values=[]
                for C0 in covs:
                    C=TT[n]@C0@TT[n].transpose(0,2,1)
                    num=(C[:,1,1]+w2*C[:,0,0])/(2*w)-.5
                    values.append(num)
                result[n,:,0]=values[0];result[n,:,1]=values[1];result[n,:,2]=w
                energy=[float(np.trapezoid(ks[ks>=.01]**3*w[ks>=.01]*val[ks>=.01]/(2*np.pi**2),np.log(ks[ks>=.01]))/EX) for val in values]
                ratios.append(energy)
                determinant=np.linalg.det(TT[n]);maxsym=max(maxsym,float(abs(determinant-1).max()))
            ratios=np.array(ratios);intervals=[]
            for end in (3.,10.,40.,80.):
                m=times<=end;intervals.append(dict(end_tau=end,maximum_conditional_phi_energy_fraction=ratios[m].max(axis=0)))
            crossing=[]
            for j in range(2):
                inds=np.flatnonzero(ratios[:,j]>=1);crossing.append(float(times[inds[0]]) if len(inds) else None)
            rows.append(dict(tag=r['tag'],intervals=intervals,first_sampled_energy_fraction_one=crossing,
                max_mode_occupation=float(result[:,:,:2].max()),minimum_occupation=float(result[:,:,:2].min()),check=check))
            profiles[r['tag']]=result.astype(np.float32);profiles[r['tag']+'_energy_ratios']=ratios
            print('NEUTRAL',rows[-1],flush=True)
    # Short independent Radau mode control on the same dynamical background.
    r=base['rows'][0];bg,tt,_=run(neutral,r,[2.,6.],[0.,3.],rtol=2e-12,suffix='_control')
    v,_=prepared(r['beta']);seed=np.array([.4,-.7])
    def fun(t,y):
        out=np.empty(21);out[:19]=rhs(t,y[:19],r['beta'],r['Gamma_over_H'],XI,ETA2)
        a=np.exp(y[15]);out[19]=y[20]/a;out[20]=-(4+mass(y[:19],r['beta']))*y[19]/a;return out
    sol=solve_ivp(fun,(0,3),np.r_[v,seed],method='Radau',rtol=2e-10,atol=2e-13)
    err=float(np.max(abs(sol.y[19:,-1]-tt[-1,0]@seed)/(1+abs(tt[-1,0]@seed))))
    np.savez_compressed(ROOT/'NEUTRAL_PROFILES.npz',**profiles)
    write('NEUTRAL_RESULTS.json',dict(rows=rows,ks=ks,seed_labels=['instantaneous canonical vacuum','first-order initial WKB seed'],
        projection='Instantaneous positive-frequency occupation on the prescribed background, with its zeroth-order zero-point energy subtracted; not a renormalized interacting stress tensor or final relic.',
        controls=dict(symplectic_determinant_error=maxsym,independent_Radau_error=err),
        scope='xi_phi=0 at matching, phi background0; no decays, scattering, backreaction or nonlinear symmetry breaking. Energy fraction reaching1 marks failure of this linear approximation, not a computed energy excess in a physical solution.'))
    verified=json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text());verified['groups']['neutral_binding_field_modes']=dict(passed=maxsym<1e-6 and err<1e-7,determinant_error=maxsym,Radau_error=err)
    verified['passed_groups']=sum(g['passed'] for g in verified['groups'].values());write('VERIFICATION_RESULTS.json',verified)
    if not verified['groups']['neutral_binding_field_modes']['passed']:raise SystemExit(1)

if __name__=='__main__':main()
