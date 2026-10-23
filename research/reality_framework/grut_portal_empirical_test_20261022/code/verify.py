"""Numerical controls and a same-author fresh replay, not independent review."""
import argparse,hashlib,importlib.metadata,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
from scipy.integrate import quad,simpson
from scipy.special import kv
from model import ROOT,PARENT,Object,bin_integrals,MH,VH,LAM,G,HD,MPLANCK
from history import ME,spectral_production,theta,DM_ENTROPY


def check_controls(cache):
 cache=Path(cache)
 provenance=json.loads((ROOT/'SOURCE_INPUTS.json').read_text())
 for path,meta in provenance['sources'].items():
  p=cache/path
  assert p.stat().st_size==meta['bytes']
  assert hashlib.sha256(p.read_bytes()).hexdigest()==meta['sha256'],path
 sys.path.insert(0,str(cache/'xenon'))
 from xe_likelihood.binwise_inference import BinwiseInference
 from xe_likelihood.spectrum import Spectrum
 det=json.loads((ROOT/'DETECTOR_RESULTS.json').read_text())
 shape=json.loads((ROOT/'SHAPE_RESULTS.json').read_text())
 pub=np.load(ROOT/'data/public_inputs.npz')
 control=BinwiseInference.from_xenon1t_sr(Spectrum.from_wimp(50))
 public_control={'accepted':float(control.sum_spectrum_in_reconstructed_bins),
  'sigma90_asym_cm2':float(control.compute_ul(asymptotic=True)*1e-45),
  'sigma90_recommended_cm2':float(control.compute_ul(asymptotic=False)*1e-45)}
 assert abs(public_control['accepted']/det['control_SI50']['accepted_events_at_1e45']-1)<1e-7
 for k in ['sigma90_asym_cm2','sigma90_recommended_cm2']:
  assert abs(public_control[k]/det['control_SI50'][k]-1)<5e-4
 obj=Object(39944);ref=next(x for x in det['cases'] if x['Q']==39944 and x['v_phi_GeV']==1.)
 observed=BinwiseInference.from_xenon1t_sr(lambda e:obj.rate(e,1.))
 obs_ul=float(observed.compute_ul(asymptotic=False))
 assert abs(obs_ul/(ref['xenon1t_kappa90_recommended']/1e-9)**2-1)<5e-4
 future=BinwiseInference.from_xenonnt_mc(lambda e:20*obj.rate(e,1.))
 toy_ul=[float(future.compute_ul(run=i)) for i in range(3)]
 assert np.allclose(toy_ul,det['reference_future_first3_ul_strength'],rtol=1e-9,atol=1e-11)
 future_sig=future.spectrum_in_reconstructed_bins
 saved=np.load(ROOT/'data/reconstructed_spectra.npz')['Q39944_v1_future_reconstructed']
 assert np.allclose(future_sig,saved,rtol=1e-8,atol=1e-10)
 # Independent boosted two-body inverse-energy moment, against the K0 result.
 # p=m*sinh(y): after the angular average, integrate 2*m*y*sinh(y).
 temp=40.;z=MH/temp
 numerical=quad(lambda y:2*MH*y*np.sinh(y)*np.exp(-z*np.cosh(y)),0,8,epsabs=1e-11)[0]
 closed=2*temp*kv(0,z)
 assert abs(numerical/closed-1)<1e-11
 # Spin-averaged ee -> X X* cross section vs scalar spectral widths.
 k=ref['kappa_anchor'];v=1.;m=.4;s=m*m;ms2=2*LAM*v*v
 be=np.sqrt(1-4*ME*ME/s);bx=np.sqrt(1-4*G*G*v*v/s)
 direct=k*k*ME*ME*G**4*v**4/(8*np.pi*(s-MH*MH)**2*(s-ms2)**2)*be*bx
 gee=ME*ME*m/(8*np.pi*VH*VH)*be**3
 gxx=G**4*v*v/(4*np.pi*m)*bx
 spectral=4*np.pi*gee*(k*v*VH)**2*gxx/((s-MH*MH)**2*(s-ms2)**2*be*be)
 assert abs(direct/spectral-1)<1e-13
 # Recoil integration conventions differ from the parent in two declared ways.
 energy=np.linspace(1,50,1001)
 rates={name:float(simpson(obj.rate(energy,1.,natural=nat,legacy=leg),x=energy))
        for name,nat,leg in [('current',True,False),('mono131_current_Helm',False,False),('natural_parent_Helm',True,True)]}
 # A very generous UV scalar spectral width, Gamma_SM=m, over 2..1000 GeV.
 # Regulate the Higgs pole with its supplied SM width; the thermal tail is
 # already exponentially tiny at 2 GeV for the stated 20 MeV maximum.
 mg=np.unique(np.r_[np.linspace(2,10,10001),np.geomspace(10,1000,1001)])
 den=(mg*mg-MH*MH)**2+(MH*.00407)**2
 gx=G**4*v*v/(4*np.pi*mg)*np.sqrt(1-4*G*G*v*v/mg**2)
 gd=HD**4*v*v/(8*np.pi*mg)*np.sqrt(1-4*(G*G+HD*HD)*v*v/mg**2)
 gh=k*k*VH*VH/(32*np.pi*mg)*((mg*mg+2*ms2)/(mg*mg-ms2))**2*np.sqrt(1-4*ms2/mg**2)
 spec=mg*(gh/den+(k*v*VH)**2*(gx+gd)/(den*(mg*mg-ms2)**2))
 n=np.linspace(0,np.log(.02/1e-5),401);temp_grid=.02*np.exp(-n)
 entropy=2*np.pi**2*10.75*temp_grid**3/45
 hubble=np.sqrt(np.pi**2*10.75/90)*temp_grid**2/MPLANCK
 power=temp_grid/np.pi**3*simpson(mg[None,:]**5*spec[None,:]*kv(2,mg[None,:]/temp_grid[:,None]),x=mg,axis=1)
 tail=float(simpson(power/(entropy*hubble),x=n))
 assert tail/DM_ENTROPY<1e-30
 photo=json.loads((ROOT/'PHOTODISINTEGRATION_RESULTS.json').read_text())
 import acropolis
 for name,expected in photo['source_sha256'].items():
  assert hashlib.sha256((Path(acropolis.__file__).parent/name).read_bytes()).hexdigest()==expected,name
 shape_data=np.loadtxt(ROOT/'data/reference_shape_comparison.csv',delimiter=',',skiprows=1)
 assert not np.any((shape_data[:,1]>0)&(shape_data[:,2]==0))
 a=photo['cases'][2]['relative_changes'];b=photo['cases'][3]['relative_changes']
 refinement={key:float((b[key]-a[key])/a[key]) for key in ['D_over_H_relative','He3_over_D_relative']}
 assert max(abs(x) for x in refinement.values())<.01
 result={'scope':'One-author numerical controls using original public implementations and distinct formula routes; not independent scientific review.',
  'official_unmodified_SI50_control':public_control,'official_reference_observed_UL_strength_at_1e9':obs_ul,
  'official_reference_future_first3_UL_strength':toy_ul,
  'boosted_decay_variance_relative_residual':float(numerical/closed-1),
  'ee_to_XX_direct_vs_spectral_relative_residual':float(direct/spectral-1),
  'recoil_convention_1_to_50keV_rates_per_tonne_year_at_1e9':rates,
  'cold_2_to_1000GeV_tail_GeV_per_entropy_with_GammaSM_le_mass':tail,
  'cold_2_to_1000GeV_tail_over_DM':tail/DM_ENTROPY,
  'ACROPOLIS_refined_relative_change_of_increment':refinement,
  'shape_fit_KL_nonnegative':bool(shape['KL_object_to_best_point_per_signal_event']>=0),
  'environment':{'python':sys.version.split()[0],**{name:importlib.metadata.version(name) for name in ['numpy','scipy','matplotlib','tqdm','acropolis','numba','llvmlite']}}}
 (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps(result,indent=2))


def numeric_leaves(value,path=''):
 if isinstance(value,dict):
  for key,item in value.items():yield from numeric_leaves(item,path+'/'+key)
 elif isinstance(value,list):
  for i,item in enumerate(value):yield from numeric_leaves(item,path+'/'+str(i))
 elif isinstance(value,(int,float)) and not isinstance(value,bool):yield path,float(value)


def replay(cache):
 with tempfile.TemporaryDirectory(prefix='grut-portal-replay-') as task:
  base=Path(task);packet=base/ROOT.name;parent=base/PARENT.name
  shutil.copytree(ROOT,packet,ignore=shutil.ignore_patterns('__pycache__','MANIFEST.json'))
  (parent/'data').mkdir(parents=True);shutil.copy2(PARENT/'RESULTS.json',parent/'RESULTS.json')
  for q in [19972,39944,79988]:shutil.copy2(PARENT/'data'/f'profile_Q{q}.csv',parent/'data'/f'profile_Q{q}.csv')
  env=dict(os.environ,OPENBLAS_NUM_THREADS='1')
  for script,args in [('detector.py',[str(Path(cache).resolve())]),('history.py',[]),('shape.py',[])]:
   subprocess.run([sys.executable,str(packet/'code'/script),*args],env=env,check=True)
  counts={}
  for name in ['DETECTOR_RESULTS.json','HISTORY_RESULTS.json','SHAPE_RESULTS.json']:
   old=dict(numeric_leaves(json.loads((ROOT/name).read_text())))
   new=dict(numeric_leaves(json.loads((packet/name).read_text())))
   assert old.keys()==new.keys()
   for key in old:assert np.isclose(old[key],new[key],rtol=3e-10,atol=1e-35),(name,key,old[key],new[key])
   counts[name]=len(old)
  for name in ['projected_toy_limits.csv','reference_cold_history.csv','reference_hot_history.csv','reference_observed_curve.csv','reference_shape_comparison.csv']:
   assert np.allclose(np.loadtxt(ROOT/'data'/name,delimiter=',',skiprows=1),np.loadtxt(packet/'data'/name,delimiter=',',skiprows=1),rtol=3e-10,atol=1e-35),name
  for name in ['reconstructed_spectra.npz','public_inputs.npz']:
   old=np.load(ROOT/'data'/name);new=np.load(packet/'data'/name)
   assert old.files==new.files
   for key in old.files:assert np.allclose(old[key],new[key],rtol=3e-10,atol=1e-35),(name,key)
  result={'same_author':True,'numeric_leaves_reproduced':counts,'CSV_tables_reproduced':5,'NPZ_archives_compared_by_arrays':2,
          'photodisintegration_replayed':False,'reason':'Separate original ACROPOLIS runs plus effective-grid refinement are archived; no duplicate full nuclear run in this replay.'}
  (ROOT/'REPLAY_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('cache');p.add_argument('--replay',action='store_true');args=p.parse_args()
 replay(args.cache) if args.replay else check_controls(args.cache)
