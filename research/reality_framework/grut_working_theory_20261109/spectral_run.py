"""Only the frozen suite invokes this physical-source calculation."""
import argparse,json
from pathlib import Path
from spectral_model import evaluate_study
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise FileExistsError('Preserve previous result')
    card=json.loads((HERE/'spectral_protocol.json').read_text())
    result=evaluate_study(card,HERE.parent/'grut_response_separation_20261109/r1_reference.py')
    result['independent_R1_check']['reference_file']='research/reality_framework/grut_response_separation_20261109/r1_reference.py'
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'numerically_valid':result['numerically_valid'],'uniform_certificate':result['primary_mathematical_uniform_certificate']}))
    raise SystemExit(0 if result['numerically_valid'] else 1)
