"""Direct nonlinear test of the original relative-shift preparation.

The unchanged Cartesian solver already supplies --impact. This additive
protocol measures actual core positions and shapes rather than interpreting
an unconverged canonical tangent norm as internal breakup.
"""
from pathlib import Path
import hashlib, json, sys, time
import run
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
pp=P/'TRANSVERSE_PROTOCOL.json'
if not pp.exists():
    cs=[]
    for tag,h in [('coarse',.75),('medium',.5),('fine',.375)]:
      cs.append(dict(id='impact_'+tag,case='pair',h=h,L=24,Lz=48,T=160,extra={'--impact':.02}))
    cs.extend([
      dict(id='impact_half_fine',case='pair',h=.375,L=24,Lz=48,T=160,extra={'--impact':.01}),
      dict(id='impact_box',case='pair',h=.5,L=30,Lz=60,T=160,extra={'--impact':.02}),
      dict(id='impact_dt',case='pair',h=.5,L=24,Lz=48,T=160,dt=.025,extra={'--impact':.02})])
    spec={'local_freeze_unix':time.time(),'public_preregistration':False,'source_sha256':sha(Path(__file__)),'solver_sha256':sha(P/'cartesian.cpp'),'parent_protocol_sha256':sha(P/'protocol.json'),'cases':cs,'posthoc_scope':'Proposed after original large-domain mode fits showed >99.6% of local geometric m1 response as core translations, with material global symmetry error. No finite pure-impact Cartesian result inspected before this freeze.','hypothesis':'Small opposite source-center displacements, with velocity directions held fixed, produce chiefly outgoing trajectory deflection. Cores retain charge and bounded intrinsic shape; centroid response should be linear in initial impact for .01/.02 and converge in space/time/box.','observables':'Core x-centroid/impact, late slope, R12 charge retention, R8 charge covariance eigenvalues. Compare unshifted controls and the original axisymmetric linear core-position fits; do not equate R36 Cartesian canonical norms with whole-domain cylindrical norms.','criteria':{'retention':'Same late tau120..160, >=90% initial Q in two nonoverlapping R12 cores, >=15% in each.','linearization':'Fine centroid response x/impact changes <2% of response magnitude when impact .02 is halved.','space':'Max late centroid-response difference medium/fine <5% of response magnitude and smaller than coarse/medium difference.','time_box':'Each max late centroid-response change <5% of response magnitude.','shape':'Report actual core covariance/charge changes vs matched reference. Bounded response over a finite interval is not a measured asymptotic Lyapunov exponent.'}}
    pp.write_text(json.dumps(spec,indent=2)+'\n')
spec=json.loads(pp.read_text())
if spec['source_sha256']!=sha(Path(__file__)) or spec['solver_sha256']!=sha(P/'cartesian.cpp'):raise RuntimeError('frozen transverse source changed')
if __name__=='__main__' and '--freeze-only' not in sys.argv:
    run.build();requested=set(sys.argv[1:])
    for c in spec['cases']:
      if requested and c['id'] not in requested:continue
      if (run.RAW/(c['id']+'.execution.json')).exists():raise RuntimeError('do not overwrite '+c['id'])
      run.execute(c)
