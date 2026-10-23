#!/usr/bin/env python3
"""Fetch pinned public inputs; optional original hadronic-table regeneration."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = 'https://files.pythonhosted.org/packages/3b/b3/528b1f2ec96ebae64d87a55d9247362834589393b9772d86d1f2486ac576/hipsofcobra-0.0.8.tar.gz'
TOY = 'xe_likelihood/data/xent_nr_20ty_projection/XENONnT_binwise_NR_LL_toyMC.json.gz'


def fetch(url, destination, digest):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        data = urllib.request.urlopen(url, timeout=60).read()
        assert hashlib.sha256(data).hexdigest() == digest, url
        destination.write_bytes(data)
    assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, destination


def main(cache, regenerate=False):
    cache = Path(cache).resolve()
    meta = json.loads((ROOT/'SOURCE_INPUTS.json').read_text())
    fetch('https://raw.githubusercontent.com/XENON1T/xenon1t_approximate_nuclear_recoil_likelihood/'
          +meta['xenon_release_commit']+'/'+TOY, cache/TOY,
          '9874ff4f0d5b2e3e0c48e497d1e685195baf7c16211645abf6e6dd921ab02885')
    archive = cache/'hipsofcobra-0.0.8.tar.gz'
    fetch(ARCHIVE, archive, meta['hipsofcobra_source_archive_sha256'])
    for name, source in meta['FORESEE_files'].items():
        fetch(source['url'], cache/name, source['sha256'])
        data = (cache/name).read_bytes()
        if source['packaging_operation'] == 'append_one_newline':
            data += b'\n'
        else:
            assert source['packaging_operation'] == 'identity'
        assert hashlib.sha256(data).hexdigest() == source['packaged_sha256']
        assert data == (ROOT/'data'/name).read_bytes()
    print('All pinned external source hashes pass.', flush=True)
    if regenerate:
        import numpy as np
        with tempfile.TemporaryDirectory(prefix='grut-hips-regenerate-') as tmp:
            tmp = Path(tmp)
            with tarfile.open(archive) as tar:
                tar.extractall(tmp, filter='data')
            sys.path.insert(0, str(tmp/'hipsofcobra-0.0.8'))
            os.chdir(tmp)
            import hipsofcobra.classes as c
            original_meta = json.loads((ROOT/'data/hipso_provenance.json').read_text())
            for path in Path(c.hipsofcobra_path).rglob('*'):
                if path.is_file() and path.name in original_meta['source_hashes']:
                    assert hashlib.sha256(path.read_bytes()).hexdigest() == original_meta['source_hashes'][path.name]
            rows=[]
            for method in ('DGL','BTPZ'):
                for particle in ('pi','K'):
                    c.rng = np.random.default_rng(meta['hipsofcobra_rng_seed'])
                    h = c.HipsofCobra(clist=[1,1,1], Pname=particle, method=method)
                    table = np.array(h.width_band_list)
                    table[:,0] = np.sqrt(table[:,0])
                    saved=np.loadtxt(ROOT/'data'/f'hipso_{particle}_{method}.csv', delimiter=',',skiprows=1)
                    assert np.allclose(table,saved,rtol=1e-12,atol=1e-25)
                    rows.append(dict(method=method,particle=particle,rows=len(table),
                                     max_absolute_difference=float(np.max(abs(table-saved)))))
                    print('Regenerated',method,particle,'agrees.',flush=True)
        (ROOT/'HADRONIC_REGENERATION_RESULTS.json').write_text(json.dumps(dict(
            scope='Original pinned hipsofcobra regenerated with the original seed; same-author source reproduction.',
            passed=True,rows=rows),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('cache',type=Path)
    parser.add_argument('--regenerate-hadronic',action='store_true')
    args=parser.parse_args()
    main(args.cache,args.regenerate_hadronic)
