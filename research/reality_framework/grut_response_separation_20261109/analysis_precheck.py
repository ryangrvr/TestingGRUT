"""Synthetic-only sanity checks of R1 analysis; never loads lattice states."""
import sys,json,copy,argparse
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
import analyze as A
from r1_reference import encode_matrix,matrix_diagnostics

times=(0.,.29,.83,1.57);freq=np.array([.7,1.3,2.1])
B=np.array([[.3,.08,-.03],[.05,.27,.04],[-.02,.06,.22],[.13,.11,.15]])
C=B@B.T;K=(B*freq)@B.T;names=['sine1','sine2','sine3','beta_left'];means=[.2,.0,.06,.17]
responses={f'{t:g}':encode_matrix((B*np.exp(-1j*freq*t))@B.T) for t in times}
primary={t:np.asarray(m)[:3,:3,:] for t,m in responses.items()}
D=matrix_diagnostics(C[:3,:3],K[:3,:3],primary)
def block(double=False):
    statep={'++':.25,'+-':.25,'-+':.25,'--':.25} if double else {'vacuum':1.,'kick':1.,'plus':.5,'minus':.5}
    dim=4 if double else 2
    return {'gram':encode_matrix(np.eye(dim)), 'probes':{
      name:{'first_matrix':encode_matrix(np.eye(dim)*means[i]),
      'second_matrix':encode_matrix(np.eye(dim)*(means[i]**2+C[i,i])),
      'states':{s:{'probability':p,'mean':means[i],'second_moment':means[i]**2+C[i,i],'variance':C[i,i]} for s,p in statep.items()}}
      for i,name in enumerate(names)}}
instant={'single':{'0':block()},'double':block(True),
  'probe_covariance':{name:{other:float(C[i,j]) for j,other in enumerate(names)} for i,name in enumerate(names)}}
ref={'means':means,'covariance':C.tolist(),'first_moment':K.tolist(),'fsum':(2*K).tolist(),
     'responses':responses,'primary_diagnostics':D}
single={s:{'probability':p,'fsum':(2*K).tolist()} for s,p in {'vacuum':1.,'kick':1.,'plus':.5,'minus':.5}.items()}
double={s:{'probability':.25,'fsum':(2*K).tolist()} for s in ('++','+-','-+','--')}
d={'probe_names':names,'means':means,'covariance':C.tolist(),'first_moment':K.tolist(),
  'responses':responses,'diagnostics':D,'refined_reference':ref,'instant_observables':instant,
  'instant_refined_reference':instant,'single_current_rules':single,'coincident_double_current_rules':double,
  'delayed_double_current_rules':double,'twice_kick_current_rule':{'double_kick':{'probability':1.,'fsum':(2*K).tolist()}},
  'instant_kick_scan':{'0.5':{'fsum':(2*K).tolist(),'maximum_link_phase':.1}}}
result=A.per_point(d,{})
checks={key:result[key] for key in ('static_covariance_relative','static_first_moment_relative',
  'max_normalized_auto_response_error','centroid_relative_error','cosine_eigen_absolute_error','decomposition_closure_error')}
checks['instant_moments']=result['instantaneous_moments']['maximum']
checks['source_span_predictor']=result['max_lattice_galerkin_residual']
# Congruence-invariant whitened relative discrepancy, including nonorthogonal mixing.
W=np.array([[1.3,.2,-.1],[.1,.8,.3],[-.2,.1,1.1]])
target=2*K[:3,:3];changed=target+np.diag([.003,-.002,.001])
e0=A.relative_whitened_error(changed,target)['operator_norm']
e1=A.relative_whitened_error(W@changed@W.T,W@target@W.T)['operator_norm']
checks['relative_current_congruence_invariance']=abs(e0-e1)
gain=copy.deepcopy(d);gain['covariance']=(1.02*C).tolist();gain['first_moment']=(1.02*K).tolist()
gain['responses']={t:encode_matrix(1.02*np.asarray(A.decode_matrix(m))) for t,m in responses.items()}
gain['diagnostics']=matrix_diagnostics(1.02*C[:3,:3],1.02*K[:3,:3],{t:np.asarray(m)[:3,:3,:] for t,m in gain['responses'].items()})
gm=A.per_point(gain,{})
checks['pure_gain_shape_removed']=gm['max_normalized_auto_response_error']
checks['pure_gain_static_recognized']=abs(gm['static_covariance_relative']-.02)
checks['both_signed_decompositions']=gm['decomposition_closure_error']
# The zero-tail Gaussian bound and normalization checker use only toy metadata.
toyref=copy.deepcopy(instant);toyref.update({'detector':{'variance':.1,'tail_variance_bound':0.},
  'probe_metadata':{name:{'variance':float(C[i,i]),'tail_variance_bound':0.,'mean':means[i]} for i,name in enumerate(names)},
  'probe_covariance_tail_bounds':{name:{other:0. for other in names} for name in names}})
checks['zero_tail_bound']=A.instant_tail_control(toyref,{'coupling':2.})['maximum']
checks['toy_normalization']=A.normalization_checks(instant)['maximum_error']
out={'scope':'Synthetic source matrices only; no physical reference or production checkpoint evaluated.',
 'checks':checks,'maximum_identity_error':max(checks.values()),'passed':max(checks.values())<1e-12}
with args.output.open('x',encoding='utf-8') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
