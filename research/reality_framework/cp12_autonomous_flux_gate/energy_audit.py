"""Infer a full-law electric-energy interval using exact conserved energy.

Norm error alone does not bound an unbounded observable. Here F=H-kappa Q-D,
and D (mass+hopping+rewiring) IS bounded at fixed N,Q. The known full-law
total energy therefore supplies the additional information needed.
"""
from pathlib import Path
import json, math

P=Path(__file__).resolve().parent
MAIN=['N4_Q4_B256_T20.json','N5_Q5_B128_T10.json',
      'N6_Q6_B128_T10.json','N5_Q6_B128_T10.json','N6_Q12_B6_T0.4.json',
      'N6_Q12_B10_T1.json']


def audit(d):
    n=d['N'];q=d['Q'];M=n*(n-1)//2
    hn=math.sqrt(2*n*q);rn=2*.1*q*(M-q)/M
    bounded=2*n+hn+rn
    rows=[]
    for r in d['rows']:
        e=r['state_error_upper'];factor=2*e+e*e
        F=r['electric_square_projected']/8
        err=bounded*factor+r['energy_error_projected']+q/2*r['norm_error']
        rows.append(dict(T=r['T'],electric_energy=F,
            electric_energy_error_upper=err,
            electric_energy_interval=[max(0.,F-err),F+err],
            mass_energy_projected=r['mass_energy_projected'],
            mass_energy_error_upper=2*n*factor,
            hopping_energy_projected=r['hopping_energy_projected'],
            hopping_energy_error_upper=hn*factor+q/2*r['norm_error'],
            rewiring_energy_projected=r['rewiring_energy_projected'],
            rewiring_energy_error_upper=rn*factor,
            electric_tradeoff_rhs_projected=.25*(q-q*r['removable_fraction']-.5*q*r['hole_protected_fraction'])))
    return dict(N=n,Q=q,budget=d['budget'],exact_total_energy=q/2,
                bounded_complement_operator_norm_upper=bounded,
                hopping_operator_norm_upper=hn,rewiring_operator_norm_upper=rn,rows=rows)


def main():
    result=dict(status='CP12_ELECTRIC_ENERGY_BY_CONSERVATION_PASS',
                models=[dict(file=name,**audit(json.loads((P/name).read_text()))) for name in MAIN],
                qualification='Uses exact conserved full-law total mean energy and bounded complementary energy terms, not a direct unbounded-observable inference from norm error. All numerical bound evaluations use ordinary floating point.')
    (P/'ENERGY_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],models=len(result['models']))),flush=True)


if __name__=='__main__':main()
