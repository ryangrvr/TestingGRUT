"""Reproduce/collect exploratory CLASS matter-power comparisons.

Use CLASS v3.4.1 commit 814ce59c380a5cccd0bed1747a384ff7bfc1f9cf, build
`make class`, then `python code/class_power.py --run --class-binary /path/class`.
The actual-spectrum hierarchy is not converged: inspect psd_audit.json and the
numerical controls. Collisionless adiabatic perturbations are additional model
assumptions and conflict with a naive claim that gas scattering is negligible.
"""
from pathlib import Path
import argparse,json,os,subprocess,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'.work'/'class_runs'
PARENT='814ce59c380a5cccd0bed1747a384ff7bfc1f9cf'
CASES={
 'cold':{'warm':False,'kmax':20},'coldtight':{'warm':False,'kmax':20,'tol':1e-7},
 'coldrk':{'warm':False,'kmax':2,'evolver':0},'coldrk20':{'warm':False,'kmax':20,'evolver':0},
 'coldhigh':{'warm':False,'kmax':1000,'evolver':0,'kdec':5},
 'warm96fast20':{'warm':True,'psd':'psd_96.dat','kmax':20,'qbins':30,'trigger':.1,'baseline':'cold'},
 'warm96q60':{'warm':True,'psd':'psd_96.dat','kmax':20,'qbins':60,'trigger':.1,'baseline':'cold'},
 'warm96rk':{'warm':True,'psd':'psd_96.dat','kmax':2,'qbins':30,'evolver':0,'baseline':'coldrk'},
 'warm96rkq60':{'warm':True,'psd':'psd_96.dat','kmax':2,'qbins':60,'evolver':0,'baseline':'coldrk'},
 'warm96rkq120':{'warm':True,'psd':'psd_96.dat','kmax':2,'qbins':120,'evolver':0,'baseline':'coldrk'},
 'warm96maxwell':{'warm':True,'psd':'psd_maxwell96.dat','kmax':20,'baseline':'cold'},
 'warm96maxwellrk':{'warm':True,'psd':'psd_maxwell96.dat','kmax':20,'qbins':30,'evolver':0,'baseline':'coldrk20'},
 'warm96maxwellhigh':{'warm':True,'psd':'psd_maxwell96.dat','kmax':1000,'evolver':0,'kdec':5,'baseline':'coldhigh'}}

def config(name):
 case=CASES[name];fs=json.loads((ROOT/'free_stream.json').read_text());f=fs['populations']['96']['mass_fraction_of_cold_plus_held'];ai=fs['background']['a_initial']
 c=json.loads((ROOT.parent/'grut_matched_candidate_20261027'/'CANDIDATE_M1.json').read_text());m=c['light_renormalized_targets']['mX_GeV']*1e9;T=m*ai*.05/(2.7255*8.617333262145e-5)
 items={'h':.674,'omega_b':.0224,'omega_cdm':.12*(1-f) if case['warm'] else .12,'N_ur':3.044,'YHe':.245,'T_cmb':2.7255,'A_s':2.1e-9,'n_s':.965,
        'recombination':'RECFAST','output':'mPk','P_k_max_1/Mpc':case['kmax'],'z_pk':0,'root':str(D/(name+'_')),'background_verbose':1,'input_verbose':1}
 if case['warm']:
  items.update({'N_ncdm':1,'use_ncdm_psd_files':1,'ncdm_psd_filenames':str(ROOT/case['psd']),'m_ncdm':m,'omega_ncdm':.12*f,'T_ncdm':T,'a_ini_over_a_today_default':1e-16})
 if 'evolver' in case:items['evolver']=case['evolver']
 if 'qbins' in case:items.update({'ncdm_quadrature_strategy':3,'ncdm_N_momentum_bins':case['qbins'],'ncdm_maximum_q':4})
 if 'trigger' in case:items['ncdm_fluid_trigger_tau_over_tau_k']=case['trigger']
 if 'tol' in case:items['tol_perturbations_integration']=case['tol']
 if 'kdec' in case:items['k_per_decade_for_pk']=case['kdec']
 return items

def run_case(name,binary):
 D.mkdir(parents=True,exist_ok=True);params=config(name);params['root']=str(D/(name+'_reproduce_'));ini=D/(name+'_reproduce.ini');ini.write_text(''.join(f'{k} = {v}\n' for k,v in params.items()))
 env=os.environ.copy();env['OMP_NUM_THREADS']='2';r=subprocess.run([str(binary),str(ini)],capture_output=True,text=True,env=env)
 (D/(name+'_reproduce.log')).write_text(r.stdout+r.stderr)
 if r.returncode:raise RuntimeError(f'{name}: CLASS status {r.returncode}; see {name}_reproduce.log')
 outputs=sorted(D.glob(name+'_reproduce_*_pk.dat'),key=lambda p:p.stat().st_mtime_ns)
 if not outputs:raise RuntimeError('CLASS did not write the requested spectrum')
 shutil.copyfile(outputs[-1],D/(name+'_00_pk.dat'))

def collect():
 runs={};arrays={};ks=(.01,.1,.3,1.,2.,5.,10.,20.)
 for name,case in CASES.items():
  f=D/(name+'_00_pk.dat')
  if not f.exists():continue
  d=np.loadtxt(f);arrays[name]=d;runs[name]={'settings':case,'points':len(d),'k_min_h_Mpc':float(d[0,0]),'k_max_h_Mpc':float(d[-1,0])}
 for name,case in CASES.items():
  if name not in arrays or not case.get('warm') or case['baseline'] not in arrays:continue
  w=arrays[name];b=arrays[case['baseline']];compare={}
  for kval in ks:
   if not (max(w[0,0],b[0,0])<=kval<=min(w[-1,0],b[-1,0])):continue
   f=lambda a:np.interp(np.log(kval),np.log(a[:,0]),np.log(a[:,1]));compare[str(kval)]=float(100*np.expm1(f(w)-f(b)))
  runs[name]['percent_difference_from_baseline']=compare
 out={'class_version':'3.4.1','class_commit':PARENT,'cosmology':{'h':.674,'omega_b':.0224,'omega_dark_total':.12,'N_ur':3.044,'recombination':'RECFAST','A_s':2.1e-9,'n_s':.965,'z':0},
      'runs':runs,'scope':'Exploratory CLASS linear matter spectra using a selected one-point velocity distribution, collisionless adiabatic ncdm present from early times, and a fixed population. Momentum-bin and solver controls are not converged for the narrow distribution; positive percent differences are numerical artifacts at this scale. M1 gas collisions and core formation/perturbations are not included. Not a prediction or likelihood.'}
 np.savez_compressed(ROOT/'class_spectra.npz',**arrays);(ROOT/'class_results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run',action='store_true');ap.add_argument('--class-binary');ap.add_argument('--cases',default='cold,warm96maxwell');a=ap.parse_args()
 if a.run:
  if not a.class_binary:raise ValueError('--class-binary required')
  for name in a.cases.split(','):
   if name not in CASES:raise ValueError(name)
   run_case(name,Path(a.class_binary).resolve())
 collect()
if __name__=='__main__':main()
