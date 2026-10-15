"""Fresh, self-contained reproduction of all new calculations and checks."""
from pathlib import Path
import argparse,contextlib,hashlib,io,json,platform,tempfile,time
import numpy as np
import scipy
import bounds,checks,prepared,hosts

P=Path(__file__).resolve().parent
DEPS=['cp12_ground_connectivity/frozen_model.py','cp12_ground_connectivity/exact_destination.py',
      'cp12_coherence_locality/three_dimensional_host.py']
INPUTS=['cp12_unfiltered_record_gate/DERIVATION.md','cp12_unfiltered_record_gate/MANIFEST.json',
        'cp12_unfiltered_record_gate/CYCLE_CAPACITY.json']


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest():
    return dict(parent_commit='36138c6140062b20c60864bce12a893c91dca8eb',
        packet_hashes={f.name:digest(f) for f in sorted(P.iterdir()) if f.is_file() and f.name not in ['MANIFEST.json','VERIFICATION.json']},
        dependency_hashes={name:digest(P.parent/name) for name in DEPS},
        input_hashes={name:digest(P.parent/name) for name in INPUTS})


def validate(m):
    for name,h in m['packet_hashes'].items():assert digest(P/name)==h,('packet',name)
    for section in ['dependency_hashes','input_hashes']:
        for name,h in m[section].items():assert digest(P.parent/name)==h,(section,name)
    assert digest(P.parent/'cp12_ground_connectivity/frozen_model.py')=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'


def compare(a,b,key=''):
    if key=='seconds':return 0.
    if isinstance(a,dict):
        assert set(a)==set(b),(key,'keys')
        return max([0.]+[compare(a[k],b[k],k) for k in a])
    if isinstance(a,list):
        assert len(a)==len(b),(key,'length')
        return max([0.]+[compare(x,y,key) for x,y in zip(a,b)])
    if isinstance(a,(int,float)) and not isinstance(a,bool):
        error=abs(float(a)-float(b));assert error<=1e-9*(1+abs(float(b))),(key,a,b)
        return error
    assert a==b,(key,a,b)
    return 0.


def run(record=False):
    start=time.monotonic();mp=P/'MANIFEST.json'
    if mp.exists() and not record:validate(json.loads(mp.read_text()))
    rows=[];maximum=0.
    with tempfile.TemporaryDirectory(prefix='cp12-extensive-verify-') as temp:
        for module,filename in [(bounds,'ANALYTIC_BOUNDS.json'),(checks,'CHECKS.json'),
                                (prepared,'CHARGED_CONTROLS.json'),(hosts,'MATCHED_HOSTS.json')]:
            oldp=module.P;module.P=Path(temp)
            try:
                with contextlib.redirect_stdout(io.StringIO()):result=module.run()
            finally:module.P=oldp
            saved=json.loads((P/filename).read_text());error=compare(result,saved)
            maximum=max(maximum,error);rows.append(dict(calculation=filename,freshly_executed=True,maximum_numeric_difference=error))
            print(json.dumps(rows[-1]),flush=True)
    if record:mp.write_text(json.dumps(manifest(),indent=2)+'\n')
    validate(json.loads(mp.read_text()))
    out=dict(status='CP12_EXTENSIVE_CIRCULATION_GATE_PASS',fresh_calculations=rows,
        maximum_saved_numeric_difference=maximum,all_new_calculations_freshly_rerun=True,
        external_cache_or_wave_snapshots_required=False,source_and_data_hashes_checked=True,
        Python=platform.python_version(),NumPy=np.__version__,SciPy=scipy.__version__,
        seconds=time.monotonic()-start,
        scope='Full new packet rerun. Exact rational identities and analytic bounds are distinguished from floating-point full-H-action, norm and generic finite-evolution checks. No extensive-state dynamics or formal interval-roundoff certificate.')
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],seconds=out['seconds'])),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true')
    run(ap.parse_args().record_manifest)
