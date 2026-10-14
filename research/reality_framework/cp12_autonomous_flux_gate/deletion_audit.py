"""Residual coherent-removal capacity after autonomous internal recording.

Exact Johnson incidence identity on every fixed matter/field fiber:
D^dag D=(M-Q+1)r-(M/Gamma)R. D is a mathematical deletion map, not a
connection-number-violating physical term added to the Hamiltonian.
"""
from pathlib import Path
import json
from energy_audit import MAIN

P=Path(__file__).resolve().parent


def audit(d):
    n=d['N'];q=d['Q'];M=n*(n-1)//2;cap=M-q+1;rows=[]
    for r in d['rows']:
        nr=q*r['removable_fraction'];R=r['rewiring_energy_projected']
        numerator=cap*nr-M/.1*R
        e=r['state_error_upper'];factor=2*e+e*e
        error=q*cap*factor;den_error=q*factor
        assert numerator>=-1e-8 and numerator<=q*cap+1e-8
        amp=[max(0.,numerator-error)/(nr+den_error),
             (numerator+error)/(nr-den_error) if nr>den_error else None]
        rows.append(dict(T=r['T'],coherent_deletion_norm_square=numerator,
            coherent_deletion_norm_square_error_upper=error,
            complete_configuration_dephased_deletion_norm_square=nr,
            coherent_to_dephased_ratio=numerator/nr if nr else None,
            coherent_to_dephased_ratio_interval=amp,
            ratio_to_removability_matched_maximum=numerator/(cap*nr) if nr else None))
    return dict(N=n,Q=q,budget=d['budget'],maximum_amplification=cap,rows=rows)


def main():
    result=dict(status='CP12_RESIDUAL_DELETION_COHERENCE_AUDIT_PASS',
                models=[dict(file=f,**audit(json.loads((P/f).read_text()))) for f in MAIN],
                scope='Exact bounded-operator diagnostic in the actual unforced state. No graph record or dephasing is applied to that state. The hypothetical dephased comparator changes energy/entropy; the deletion norm is not a finite-time matter signal or a causal-cone test.')
    (P/'DELETION_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],models=len(result['models']))),flush=True)


if __name__=='__main__':main()
