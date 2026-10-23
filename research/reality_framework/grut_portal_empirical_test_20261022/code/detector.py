"""Recast official observed and projected likelihood arrays."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy.stats import chi2
from model import ROOT,Object,bin_integrals


def curve_grid(sig,ys,likelihood,n=1000):
 strength=np.linspace(0,ys[-1]/max(sig),n)
 data=likelihood if likelihood.ndim==3 else likelihood[None,:,:]
 out=np.zeros((len(data),n))
 for i,s in enumerate(sig):
  pos=s*strength
  j=np.clip(np.searchsorted(ys,pos,side='right')-1,0,len(ys)-2)
  t=(pos-ys[j])/(ys[j+1]-ys[j])
  out+=data[:,i,j]*(1-t)+data[:,i,j+1]*t
 out-=np.min(out,axis=1)[:,None]
 return strength,out


def limits(strength,ll,sigsum,threshold=None):
 target=chi2(1).isf(.1) if threshold is None else np.interp(strength*sigsum,*threshold,left=chi2(1).isf(.1),right=chi2(1).isf(.1))
 diff=ll-target
 out=[]
 for row,cur in zip(diff,ll):
  j0=int(np.argmin(cur))
  eligible=np.flatnonzero(row[j0:]>=0)
  if not len(eligible):raise ValueError('Limit outside public likelihood range')
  j=j0+eligible[0]
  if j==0:raise ValueError('Limit at lower grid edge')
  # Interpolation exactly solves the piecewise-linear likelihood used here.
  out.append(brentq(lambda mu:np.interp(mu,strength,row),strength[j-1],strength[j],xtol=1e-13))
 return np.array(out)


def run(cache):
 cache=Path(cache)
 pub=np.load(ROOT/'data/public_inputs.npz')
 futurefile=cache/'xenon/xe_likelihood/data/xent_nr_20ty_projection/XENONnT_binwise_NR_LL_toyMC.json.gz'
 provenance=json.loads((ROOT/'SOURCE_INPUTS.json').read_text())
 expected=provenance['sources'][str(futurefile.relative_to(cache))]['sha256']
 assert hashlib.sha256(futurefile.read_bytes()).hexdigest()==expected
 with gzip.open(futurefile,'rt') as f:raw=json.load(f)
 data=raw['Z']['Data'];toy=np.array([data[str(i)] for i in range(len(data))]);del raw,data
 # The public toy signal axis uses the same signal range as its JSON header.
 with gzip.open(futurefile,'rt') as f:
  header=f.read(10000)
 # Header is before the large Z field; parse its two axis objects directly.
 decoder=json.JSONDecoder();start=header.index('"Yedges"');axis=decoder.raw_decode(header[header.index('{',start):])[0]
 ye=np.round(np.linspace(axis['low_edge'],axis['high_edge'],axis['N']),6);ty=(ye[1:]+ye[:-1])/2
 oy=(pub['obs_ll_yedges'][1:]+pub['obs_ll_yedges'][:-1])/2
 threshold=(pub['threshold_x'],pub['threshold_y'])
 result={'classification':'Official XENON1T approximate likelihood recast and named historical XENONnT 20 ty projection; LZ is an efficiency-reweighting estimate only.',
         'kappa_reference':1e-9,'halo_fraction':1.,'projection_toys':len(toy),'cases':[],'checks':{}}
 print('Loaded',len(toy),'projection toys',toy.shape,flush=True)
 saved={};toycols=[];headers=['toy']
 obj0=Object(39944)
 control=interp1d(pub['control_SI50_energy'],pub['control_SI50_rate'],bounds_error=False,fill_value=0.)
 cbin=bin_integrals(control,pub['obs_response_xedges'],16)
 csig=cbin@pub['obs_response_values']
 st,ll=curve_grid(csig,oy,pub['obs_ll_values'])
 result['control_SI50']={'accepted_events_at_1e45':float(sum(csig)),
   'sigma90_asym_cm2':float(limits(st,ll,sum(csig))[0]*1e-45),
   'sigma90_recommended_cm2':float(limits(st,ll,sum(csig),threshold)[0]*1e-45)}
 for q in (19972,39944,79988):
  obj=Object(q)
  for v in (.5,.75,1.,1.25):
   fun=lambda e:obj.rate(e,v)
   ob=bin_integrals(fun,pub['obs_response_xedges']);fb=20*bin_integrals(fun,pub['future_response_xedges'])
   osig=ob@pub['obs_response_values'];fsig=fb@pub['future_response_values']
   ost,oll=curve_grid(osig,oy,pub['obs_ll_values'])
   observed=limits(ost,oll,sum(osig),threshold)[0]
   refined_st,refined_ll=curve_grid(osig,oy,pub['obs_ll_values'],4000)
   observed_refined=limits(refined_st,refined_ll,sum(osig),threshold)[0]
   fst,fll=curve_grid(fsig,ty,toy)
   projected=limits(fst,fll,sum(fsig))
   kvals=1e-9*np.sqrt(projected);pcts=np.quantile(kvals,[.1586552539,.5,.8413447461])
   # LZ's nominal WIMP limit is mapped by equal efficiency-weighted counts.
   ee=np.linspace(0.001,100,4001);eff=np.interp(ee,pub['lz_efficiency_x'],pub['lz_efficiency_y0'],left=0,right=0)
   extended=simpson(fun(ee)*eff,x=ee)
   point=simpson(obj.rate(ee,v,point=True)*eff,x=ee)
   sx=float(np.exp(np.interp(np.log(obj.mass(v)),np.log(pub['lz_si_x']),np.log(pub['lz_si_y0']))))
   sk=obj.sigma_n(v,1e-9);ratio=extended/point
   klz=1e-9*np.sqrt(sx/(sk*ratio))
   # Grid controls are numerical, not errors on the approximate statistical model.
   alt=bin_integrals(fun,pub['obs_response_xedges'],16)
   row={'Q':q,'v_phi_GeV':v,'mass_GeV':obj.mass(v),'sigma_n_at_reference_cm2':sk,
        'xenon1t_accepted_at_reference':float(sum(osig)),
        'xenon1t_kappa90_recommended':float(1e-9*np.sqrt(observed_refined)),
        'xenon1t_strength_limit_grid_relative':float(observed/observed_refined-1),
        'bin_quadrature_order16_relative_max':float(max(abs(alt-ob))/max(ob)),
        'xenonnt_projected_kappa90_percentiles':pcts.tolist(),
        'xenonnt_accepted_at_reference':float(sum(fsig)),
        'lz2024_point_sigma90_cm2':sx,'lz2024_extended_to_point_efficiency_ratio':float(ratio),
        'lz2024_count_reweighted_kappa90_ESTIMATE':float(klz),
        'nominal_projection_to_lz_coupling_ratio':float(pcts[1]/klz),
        'kappa_anchor':float(np.sqrt(pcts[1]*min(klz,1e-9*np.sqrt(observed_refined))))}
   result['cases'].append(row)
   name=f'Q{q}_v{v:g}'
   headers.append(name);toycols.append(kvals)
   saved[name+'_obs_true']=ob;saved[name+'_obs_reconstructed']=osig
   saved[name+'_future_true']=fb;saved[name+'_future_reconstructed']=fsig
   if q==39944 and v==1:
    np.savetxt(ROOT/'data/reference_observed_curve.csv',np.column_stack([ost,oll[0],np.interp(ost*sum(osig),*threshold)]),delimiter=',',header='strength,delta_minus2logL,threshold90',comments='')
    result['reference_future_first3_ul_strength']=projected[:3].tolist()
   print(name,'observed',row['xenon1t_kappa90_recommended'],'projected',pcts[1],'LZ estimate',klz,flush=True)
 np.savez_compressed(ROOT/'data/reconstructed_spectra.npz',**saved)
 np.savetxt(ROOT/'data/projected_toy_limits.csv',np.column_stack([np.arange(len(toy)),*toycols]),delimiter=',',header=','.join(headers),comments='')
 (ROOT/'DETECTOR_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('cache');run(p.parse_args().cache)
