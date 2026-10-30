"""Optional CLASS collisionless control using the post-transport smooth shape.

This extrapolates the smooth relic back before its modeled formation and omits
self-scattering. It cannot be the collisional M1 prediction.
"""
from pathlib import Path
import json,os,subprocess,argparse,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent;W=ROOT/'.work/class';W.mkdir(parents=True,exist_ok=True)
def config(tag,warm,kmax,nq,tol):
 fs=json.loads((RF/'grut_m1_synthesis_20261029/free_stream.json').read_text());f=fs['populations']['96']['mass_fraction_of_cold_plus_held'];ai=fs['background']['a_initial'];p=json.loads((RF/'grut_matched_candidate_20261027/CANDIDATE_M1.json').read_text());mass=p['light_renormalized_targets']['mX_GeV']*1e9
 d={'h':.674,'omega_b':.0224,'omega_cdm':.12*(1-f) if warm else .12,'N_ur':3.044,'YHe':.245,'T_cmb':2.7255,'A_s':2.1e-9,'n_s':.965,'recombination':'RECFAST','output':'mPk','P_k_max_1/Mpc':kmax,'z_pk':0,'root':str(W/(tag+'_')),'k_per_decade_for_pk':2,'tol_perturbations_integration':tol,'input_verbose':1,'perturbations_verbose':1}
 if warm:d.update({'N_ncdm':1,'use_ncdm_psd_files':1,'ncdm_psd_filenames':str(ROOT/'elastic_s96_n512_a16_psd.dat'),'m_ncdm':mass,'omega_ncdm':.12*f,'T_ncdm':mass*ai*.05/(2.7255*8.617333262145e-5),'a_ini_over_a_today_default':1e-16,'ncdm_quadrature_strategy':3,'ncdm_N_momentum_bins':nq,'ncdm_maximum_q':8})
 return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--binary',required=True);p.add_argument('--tag',default='coldmicro');p.add_argument('--warm',action='store_true');p.add_argument('--kmax',type=float,default=1e5);p.add_argument('--nq',type=int,default=40);p.add_argument('--tol',type=float,default=1e-5);p.add_argument('--evolver',type=int,default=1);p.add_argument('--timeout',type=float,default=240);a=p.parse_args();d=config(a.tag,a.warm,a.kmax,a.nq,a.tol);d['evolver']=a.evolver;ini=W/(a.tag+'.ini');ini.write_text(''.join(f'{k} = {v}\n' for k,v in d.items()));env=os.environ.copy();env['OMP_NUM_THREADS']='1';tic=time.monotonic()
 with (W/(a.tag+'.log')).open('w') as f:
  try:r=subprocess.run([str(Path(a.binary).resolve()),str(ini)],stdout=f,stderr=subprocess.STDOUT,env=env,timeout=a.timeout);status=r.returncode
  except subprocess.TimeoutExpired:status='timeout'
 out={'settings':vars(a),'status':status,'seconds':time.monotonic()-tic,'parameters':d,'scope':__doc__};(ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
 if status==0:
  files=sorted(W.glob(a.tag+'_*_pk.dat'),key=lambda p:p.stat().st_mtime_ns)
  if files:np.savez_compressed(ROOT/(a.tag+'.npz'),spectrum=np.loadtxt(files[-1]))
