"""Fetch pinned public inputs into an explicitly chosen external cache."""
import argparse,hashlib,json,subprocess
from pathlib import Path
from model import ROOT


def run(cache):
 cache=Path(cache).resolve();cache.mkdir(parents=True,exist_ok=True)
 meta=json.loads((ROOT/'SOURCE_INPUTS.json').read_text());xenon=cache/'xenon'
 if not xenon.exists():
  subprocess.run(['git','clone','--depth','1','--branch','v1.0.0','--filter=blob:none','--sparse',
   'https://github.com/XENON1T/xenon1t_approximate_nuclear_recoil_likelihood.git',str(xenon)],check=True)
  subprocess.run(['git','-C',str(xenon),'sparse-checkout','set','xe_likelihood'],check=True)
 commit=subprocess.check_output(['git','-C',str(xenon),'rev-parse','HEAD'],text=True).strip()
 assert commit==meta['xenon_commit'],('Unexpected source commit',commit)
 urls={
  'lz_si.json':'https://www.hepdata.net/record/data/155182/1775562/2',
  'lz_efficiency.json':'https://www.hepdata.net/record/data/155182/1775565/2'}
 for name,url in urls.items():
  path=cache/name
  if not path.exists():subprocess.run(['curl','--fail','--location','--retry','2',url,'--output',str(path)],check=True)
 for name,item in meta['sources'].items():
  path=cache/name
  assert path.stat().st_size==item['bytes'],name
  assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256'],name
 print('All nine pinned public inputs match byte counts and SHA-256 hashes.')


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('cache');run(p.parse_args().cache)
