"""Extract exact public arrays; keep the large toy ensemble in a pinned cache."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def hist(p):
    d = json.loads(p.read_text())
    x = d['Xedges']; y = d['Yedges']
    xe = np.round(np.linspace(x['low_edge'], x['high_edge'], x['N']), 6)
    ye = np.round(np.linspace(y['low_edge'], y['high_edge'], y['N']), 6)
    return xe, ye, np.asarray(d['Z']['Data']['0'])


def run(cache):
    cache = Path(cache)
    xe = cache/'xenon/xe_likelihood/data'
    obs = xe/'xe1t_nr_1ty_result'
    future = xe/'xent_nr_20ty_projection'
    paths = [obs/'XENON1T_binwise_NR_MM.json', obs/'XENON1T_binwise_NR_LL.json',
             obs/'threshold_using_lines_90pct.json', future/'XENONnT_binwise_NR_MM.json',
             future/'XENONnT_binwise_NR_LL_toyMC.json.gz',
             xe/'ancillary_data/dRdE_wimp_SI_50GeV_10-45.csv',
             xe/'ancillary_data/xenon1t_1ty.csv',
             cache/'lz_si.json', cache/'lz_efficiency.json']
    arrays = {}
    for prefix, file in [('obs_response', paths[0]), ('obs_ll', paths[1]), ('future_response', paths[3])]:
        x,y,z = hist(file)
        arrays[prefix+'_xedges'] = x
        arrays[prefix+'_yedges'] = y
        arrays[prefix+'_values'] = z
    threshold = json.loads(paths[2].read_text())
    arrays['threshold_x'] = np.array(threshold['X']['Data'])
    arrays['threshold_y'] = np.array(threshold['Y']['Data'])
    control = np.loadtxt(paths[5], delimiter=',')
    arrays['control_SI50_energy'] = control[:,0]
    arrays['control_SI50_rate'] = control[:,1]
    for name in ('lz_si','lz_efficiency'):
        d = json.loads((cache/(name+'.json')).read_text())
        arrays[name+'_x'] = np.array([float(v['x'][0]['value']) for v in d['values']])
        for j in range(len(d['values'][0]['y'])):
            arrays[name+'_y'+str(j)] = np.array([float(v['y'][j]['value']) for v in d['values']])
    (ROOT/'data').mkdir(exist_ok=True)
    np.savez_compressed(ROOT/'data/public_inputs.npz', **arrays)
    provenance = {'xenon_commit': '8254684105193d4a84cb8b5752f5e450d5ab70b3',
                  'sources': {str(p.relative_to(cache)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in paths},
                  'array_metadata': {k: {'shape': list(v.shape), 'dtype': str(v.dtype)} for k,v in arrays.items()},
                  'large_external_input': 'xenon/xe_likelihood/data/xent_nr_20ty_projection/XENONnT_binwise_NR_LL_toyMC.json.gz'}
    (ROOT/'SOURCE_INPUTS.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps({'arrays':len(arrays),'archive_bytes':(ROOT/'data/public_inputs.npz').stat().st_size}))


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('cache');run(p.parse_args().cache)
