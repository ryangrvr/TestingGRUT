"""Assess frozen criteria without converting failed science checks to passes."""
from pathlib import Path
import hashlib, json
import numpy as np

P=Path(__file__).resolve().parent;R=P/'raw'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return np.genfromtxt(R/(name+'.csv'),delimiter=',',names=True,ndmin=1)
def core(a,i):
    q=a[f'c{i}_Q8'];cov=np.empty((len(a),3,3))
    for x,y,z in [(0,0,'xx'),(0,1,'xy'),(0,2,'xz'),(1,1,'yy'),(1,2,'yz'),(2,2,'zz')]:cov[:,x,y]=cov[:,y,x]=a[f'c{i}_{z}']/q
    ev=np.linalg.eigvalsh(cov)
    return {'Q12':a[f'c{i}_Q12'],'EQ12':a[f'c{i}_E12']/a[f'c{i}_Q12'],'anisotropy_R8':np.sqrt(ev[:,-1]/np.maximum(ev[:,0],1e-30))}
def plain(x):
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    return x
def assess(c):
    a=read(c['id']);late=a[a['t']>=120];Q0=a[0]['Q'];i=int(np.argmax([a[-1]['c0_Q12'],a[-1]['c1_Q12']]))
    out={'id':c['id'],'case':c,'N':round(2*c['L']/c['h']),'Nz':round(2*c.get('Lz',c['L'])/c['h']),'initial_Q':Q0,'initial_E':a[0]['E'],'initial_EQ':a[0]['E']/abs(Q0),'last_t':a[-1]['t'],'max_relative_Q_drift':max(abs(a['Q']/Q0-1)),'max_relative_E_drift':max(abs(a['E']/a[0]['E']-1)),'endpoint_Qnegative':a[-1]['Qnegative'],'endpoint_maxF2':a[-1]['maxF2'],'endpoint_minS':a[-1]['minS']}
    out['conservation_pass']=bool(out['max_relative_Q_drift']<1e-9 and out['max_relative_E_drift']<1e-3)
    if c['case'] in ['single','pair']:
      out['late_core_Q12_means']=[np.mean(late[f'c{j}_Q12']) for j in [0,1]]
      out['late_core_EQ12_ranges']=[[min(late[f'c{j}_E12']/late[f'c{j}_Q12']),max(late[f'c{j}_E12']/late[f'c{j}_Q12'])] for j in [0,1]]
      out['late_core_anisotropy_ranges']=[[min(core(late,j)['anisotropy_R8']),max(core(late,j)['anisotropy_R8'])] for j in [0,1]]
      out['endpoint_centers']=[[a[-1][f'c{j}_{x}'] for x in ['x','y','z']] for j in [0,1]]
      if c['case']=='single':
        out['minimum_late_Q12_fraction']=min(late[f'c{i}_Q12'])/Q0
        out['retention_pass']=bool(out['minimum_late_Q12_fraction']>=.9 and max(late[f'c{i}_E12']/late[f'c{i}_Q12'])<1)
      else:
        sep=np.sqrt(sum((late[f'c0_{x}']-late[f'c1_{x}'])**2 for x in ['x','y','z']))
        out['minimum_late_Q12_fraction'] =min(late['c0_Q12']+late['c1_Q12'])/Q0
        out['minimum_late_separation']=min(sep)
        out['retention_pass']=bool(out['minimum_late_Q12_fraction']>=.9 and min(late['c0_Q12'])/Q0>=.15 and min(late['c1_Q12'])/Q0>=.15 and min(sep)>24)
    comp=np.genfromtxt(R/(c['id']+'.csv.components.csv'),delimiter=',',names=True,ndmin=1)
    out['components']={str(f):{'count':len(comp[np.isclose(comp['threshold'],f)]),'largest_Q':sorted(comp['Q'][np.isclose(comp['threshold'],f)].tolist(),reverse=True)[:12],'selected_total_Q':sum(comp['Q'][np.isclose(comp['threshold'],f)]),'selected_fraction_of_total_Q':sum(comp['Q'][np.isclose(comp['threshold'],f)])/Q0} for f in [.02,.05,.1]}
    if c['case']=='single':
      selected=out['components']['0.02']['largest_Q'];out['single_dominant_component_pass']=bool(len(selected)>0 and selected[0]>.8*Q0 and (len(selected)==1 or selected[1]<.05*Q0))
      out['retention_pass']=bool(out['retention_pass'] and out['single_dominant_component_pass'])
    if 'tan0_norm' in a.dtype.names:
      out['endpoint_tangent_amplification']=[a[-1][f'tan{j}_norm']/a[0][f'tan{j}_norm'] for j in range(3 if c.get('sym')=='all' else 2)]
      out['endpoint_symmetry_relative_residual']=[a[-1][f'res{j}_norm']/a[-1][f'gen{j}_norm'] for j in [0,1]]
      out['symmetry_trace']={str(t):[a[k][f'res{j}_norm']/a[k][f'gen{j}_norm'] for j in [0,1]] for t in [20,40,60,80,120,160] for k in [np.argmin(abs(a['t']-t))]}
      if c.get('sym')=='all':
        ortho=[]
        for r in a:
          G=np.array([r[f'gram{j}'] for j in range(9)]).reshape(3,3);v=np.array([r[f'dot{j}'] for j in range(3)])
          ortho.append(np.sqrt(max(0,r['tan2_norm']**2-v@np.linalg.pinv(G,rcond=1e-12)@v))/a[0]['tan2_norm'])
        out['relative_shift_after_global_translation_rotation_phase_projection']=ortho[-1]
        out['projection_note']='Canonical norm, R36, global rigid generators only; physical separate-core motion remains. This is not a Lyapunov exponent.'
    return out
def compare(a,b):
    x=read(a);y=read(b);common=np.intersect1d(x['t'],y['t']);common=common[common>=120]
    ix=np.searchsorted(x['t'],common);iy=np.searchsorted(y['t'],common)
    qs=[max(abs(x[f'c{j}_Q12'][ix]-y[f'c{j}_Q12'][iy]))/x[0]['Q'] for j in [0,1]]
    out={'cases':[a,b],'maximum_late_Q12_difference_over_initial_Q':qs}
    if 'tan2_norm' in x.dtype.names and 'tan2_norm' in y.dtype.names:
      gx=x['tan2_norm'][ix]/x[0]['tan2_norm'];gy=y['tan2_norm'][iy]/y[0]['tan2_norm']
      out['maximum_late_relative_shift_amplification_relative_difference']=max(abs(gx-gy)/np.maximum(abs(gx),1e-30))
    return out
def transverse(spec):
    cases={c['id']:c for c in spec['cases']};out={};responses={}
    for name,c in cases.items():
      a=read(name);delta=c['extra']['--impact'];m=a['t']>=120;t=a['t'][m];v=np.array([a[f'c{i}_x'][m]/delta for i in [0,1]])
      responses[name]=v
      ref='pair_reference_'+('box' if name=='impact_box' else 'dt' if name=='impact_dt' else 'coarse' if name=='impact_coarse' else 'medium' if name=='impact_medium' else 'fine')
      b=read(ref);ib=np.searchsorted(b['t'],t);shape=[]
      for i in [0,1]:
        aa=core(a[m],i)['anisotropy_R8'];bb=core(b[ib],i)['anisotropy_R8'];shape.append(float(max(abs(aa-bb))))
      fits=[]
      for i in [0,1]:
        slope,intercept=np.polyfit(t,v[i],1);err=max(abs(v[i]-(slope*t+intercept)))
        fits.append({'slope_per_tau':slope,'intercept':intercept,'maximum_fit_error':err,'maximum_fit_error_over_response':err/max(abs(v[i]))})
      out[name]={'endpoint_x_over_impact':v[:,-1],'late_fits':fits,'maximum_core_anisotropy_change_vs_reference':shape,'initial_Q_change_vs_reference':a[0]['Q']/b[0]['Q']-1,'initial_E_change_vs_reference':a[0]['E']/b[0]['E']-1}
    def diff(a,b):
      x=responses[a];y=responses[b]
      return np.max(abs(x-y),axis=1)/np.maximum(np.max(abs(y),axis=1),1e-30)
    cm=diff('impact_coarse','impact_medium');mf=diff('impact_medium','impact_fine')
    half=diff('impact_half_fine','impact_fine');box=diff('impact_box','impact_medium');dt=diff('impact_dt','impact_medium')
    checks={'space':bool(max(mf)<.05 and max(mf)<max(cm)),'linearization':bool(max(half)<.02),'box':bool(max(box)<.05),'time':bool(max(dt)<.05)}
    return {'cases':out,'relative_response_differences':{'coarse_medium':cm,'medium_fine':mf,'half_seed':half,'box':box,'dt':dt},'checks':checks,'sign_note':'Cartesian --impact=.02 shifts left center -delta, right +delta; the original cylindrical relative-shift tangent has the opposite sign. Compare magnitudes or reverse the latter.','inference':'A converged centroid response and retained core charge/shape support finite scattering deflection over this interval; no asymptotic internal-growth exponent inferred.'}
def main():
    spec=json.loads((P/'protocol.json').read_text())
    supplement=json.loads((P/'TRANSVERSE_PROTOCOL.json').read_text())
    cases=spec['cases']+supplement['cases']
    missing=[c['id'] for c in cases if not (R/(c['id']+'.execution.json')).exists()]
    if missing:raise RuntimeError('campaign incomplete: '+', '.join(missing))
    integrity={name:sha(P/name)==digest for name,digest in spec['files_sha256'].items()}
    integrity['transverse_source']=sha(P/'run_transverse.py')==supplement['source_sha256']
    integrity['transverse_solver']=sha(P/'cartesian.cpp')==supplement['solver_sha256']
    results={c['id']:assess(c) for c in cases}
    grid={}
    for obj in ['single_reference','single_generic','pair_reference','pair_generic']:
      cm=compare(obj+'_coarse',obj+'_medium');mf=compare(obj+'_medium',obj+'_fine')
      d0=max(cm['maximum_late_Q12_difference_over_initial_Q']);d1=max(mf['maximum_late_Q12_difference_over_initial_Q'])
      out={'coarse_medium':cm,'medium_fine':mf,'charge_convergence_pass':bool(d1<d0 and d1<.02)}
      if obj=='pair_reference':
        n0=cm['maximum_late_relative_shift_amplification_relative_difference'];n1=mf['maximum_late_relative_shift_amplification_relative_difference']
        out['relative_shift_norm_convergence_pass']=bool(n1<n0 and n1<.05)
      grid[obj]=out
    controls={}
    for name,a,b in [('single_box','single_generic_medium','single_generic_box'),('single_dt','single_generic_medium','single_generic_dt'),('pair_box','pair_reference_medium','pair_reference_box'),('pair_dt','pair_reference_medium','pair_reference_dt')]:
      d=compare(a,b);d['charge_control_pass']=max(d['maximum_late_Q12_difference_over_initial_Q'])<.01
      if 'maximum_late_relative_shift_amplification_relative_difference' in d:d['relative_shift_control_pass']=d['maximum_late_relative_shift_amplification_relative_difference']<.05
      controls[name]=d
    orientation={x:compare('pair_reference_'+x,'pair_rotated_'+x) for x in ['medium','fine']}
    orientation['decreasing_charge_difference']=max(orientation['fine']['maximum_late_Q12_difference_over_initial_Q'])<max(orientation['medium']['maximum_late_Q12_difference_over_initial_Q'])
    old=[]
    for p in sorted(R.glob('large_modes_h*.json')):
      d=json.loads(p.read_text());r=d['trace'][-1]
      old.append({'h':d['cfg']['h'],'last_t':r['time'],'endpoint_amplification':np.array(r['norm'])/d['trace'][0]['norm'],'endpoint_symmetry_residual':r['relative_symmetry_residual'],'late_core_shape':r['late_core_shape'],'wall_seconds':d['wall_seconds']})
    geometric=json.loads((P/'GEOMETRIC_COMPARISON.json').read_text());empirical=json.loads((P/'EMPIRICAL_COMPARISON.json').read_text())
    out={'integrity':integrity,'case_count':len(results),'cases':results,'grid_comparisons':grid,'time_boundary_controls':controls,'whole_event_rotation_control':orientation,'large_domain_original_m1':old,'direct_transverse':transverse(supplement),'spectral':json.loads((P/'SPECTRAL_RESULTS.json').read_text()),'empirical_geometric':geometric,'empirical_cold_matter':empirical,'claim_boundaries':['Finite generic-perturbation family, tau160, three-field D=0 flat classical sector, not universal survival or a quantum lifetime.','Numerical fixed-Q isolated-branch Hessian plus analytic angular ordering, not a rigorous global minimizer theorem.','Symmetry residuals must converge before using a raw canonical m1 norm as internal-instability evidence.','One conditional cosmological preparation is rejected; no viable empirically validated GRUT point established.']}
    if not all(integrity.values()):raise RuntimeError('frozen source integrity failure')
    (P/'SUMMARY.json').write_text(json.dumps(plain(out),indent=2,allow_nan=False)+'\n')
    print(json.dumps(plain({'cases':len(results),'integrity':integrity,'all_conservation_pass':all(r['conservation_pass'] for r in results.values()),'finite_objects_retention':{k:r['retention_pass'] for k,r in results.items() if 'retention_pass'in r},'grid':grid,'controls':controls,'rotation':orientation}),indent=2))
if __name__=='__main__':main()
