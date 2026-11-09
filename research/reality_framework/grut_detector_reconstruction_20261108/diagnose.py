"""Exploratory decomposition of the already-known B1 geometry miss only."""
import argparse,json,time
from pathlib import Path
from model import HERE,references,distances,signed_difference,half_l1,algebra_checks

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    card=json.loads((HERE/'protocol_b2_01.json').read_text())
    card['protocols']=card['protocols'][:4]
    start=time.monotonic()
    low=references(card,20,card['modes'],diagnostic=True)
    high=references(card,20,card['refined_modes'],diagnostic=True)
    names=['target','interior_only','edge_only','B1'];p={name:low[name]['histories'] for name in names}
    contrasts={
      'interior_with_true_edges':signed_difference(p['interior_only'],p['target']),
      'edge_after_interior':signed_difference(p['B1'],p['interior_only']),
      'edge_with_smooth_interior':signed_difference(p['edge_only'],p['target']),
      'interior_after_edges':signed_difference(p['B1'],p['edge_only'])}
    result={'label':'exploratory_B1_geometry_diagnosis_only','references':low,
      'contrasts':contrasts,'half_l1':{k:half_l1(v) for k,v in contrasts.items()},
      'mode_doubling_tv':{k:distances(low[k]['histories'],high[k]['histories']) for k in names},
      'algebra_checks':algebra_checks(card),'seconds':time.monotonic()-start,
      'scope':'No new lattice output, B2 physical prediction, or confirmation timing evaluated.'}
    with args.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps({k:result[k] for k in ('half_l1','algebra_checks','seconds')},indent=2))
if __name__=='__main__':main()
