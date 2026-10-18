"""Read-only exact eligibility arithmetic. No physical evolution or search."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json

J=F(1,3840);jp=F(1,64);jm=F(1,96)
quality=F(1,100);pstar=F(1,10)
N0=2**28;L0=((N0//2+23)//24)*24
formal_L=2016;formal_T=F(formal_L)/(6*J);T0=F(L0)/(6*J)
e2=4*quality*quality*pstar
access2=e2/(4*jp*jp)
current=F(21,10)+4*(jp+jm)+2*J
p0=F(514089,4000000);q0=F(434281,1000000)


def ratio(x):
    x=F(x)
    return {'numerator':str(x.numerator),'denominator':str(x.denominator)}


with localcontext() as ctx:
    ctx.prec=60
    dec=lambda x:Decimal(x.numerator)/Decimal(x.denominator)
    result={
      'role':'pre-execution feasibility arithmetic; no new physical campaign',
      'quality_amplitude_squared':ratio(e2),
      'quality_amplitude_display':str(dec(e2).sqrt()),
      'positive_access_integral_squared_limit':ratio(access2),
      'positive_access_integral_display':str(dec(access2).sqrt()),
      'formal_fixed_path_minimum_L':str(formal_L),'formal_T':ratio(formal_T),
      'positive_crossing_majorant_formal':ratio(2*jp*formal_T),
      'uniform_population_sufficient_formal':ratio(access2/formal_T**2),
      'uniform_population_display_formal':str(dec(access2/formal_T**2)),
      'inherited_minimum_N':str(N0),'smallest_eligible_whole_word_L':str(L0),
      'whole_word_T':ratio(T0),'positive_crossing_majorant_whole_word':ratio(2*jp*T0),
      'uniform_population_sufficient_whole_word':ratio(access2/T0**2),
      'uniform_population_display_whole_word':str(dec(access2/T0**2)),
      'current_norm_majorant':ratio(current),
      'integrated_current_majorant_formal':ratio(current*formal_T),
      'integrated_current_majorant_whole_word':ratio(current*T0),
      'output_probability_margin':ratio(p0-F(1,10)),
      'retention_probability_margin':ratio(q0-F(1,4)),
      'success_operation_quality_budget':ratio(p0/F(101)),
      'receiver_modes_at_minimum_N':str(64*N0*N0),
      'dense_one_boson_complex128_bytes':str(16*64*N0*N0),
      'authorized_memory_bytes':str(2**30),
      'dense_memory_ratio':str((16*64*N0*N0)//(2**30)),
      'physical_trajectories':0,'research_route_selected':None,
      'physical_task_verdict':'unresolved','decision':'specialist_review_stop'}
print(json.dumps(result,indent=2,sort_keys=True))
