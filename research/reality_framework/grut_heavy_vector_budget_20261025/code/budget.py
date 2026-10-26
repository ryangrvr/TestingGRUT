"""Vector-only and disjoint scalar energy ledgers, with finite-Q scope."""
import json,math
from background import ROOT,write


def main():
    b=json.loads((ROOT/'BACKGROUND_RESULTS.json').read_text())
    v=json.loads((ROOT/'VECTOR_RESULTS.json').read_text())
    p,l=b['parameters'],b['late'];H=p['H_GeV'];ex=l['E_X'];nx=abs(l['N_X'])
    mx=.09240896373753621
    rows=[]
    for name,e in [('specified',1e-7),('sensitivity',1e-6)]:
        massA=6*e*p['h_GeV'];width=e*e*massA/(48*math.pi)*(1-4*mx*mx/massA**2)**1.5
        decay_a=math.sqrt(math.pi*H/(2*width))
        deposit=massA/H*decay_a/ex
        spectrum=v['spectra']['actual_fine' if name=='specified' else 'e_sensitivity']
        number=spectrum['N_A_band'];R=number*deposit
        massL=math.sqrt(p['M_Y_GeV']**2-math.sqrt(2)*p['mu_GeV']*p['h_GeV'])
        widthL=3*p['kappa']**2*massL/(256*math.pi**3)
        RY=(l['N_YH']+l['N_YL'])*massL/H*math.sqrt(math.pi*H/(2*widthL))/ex
        RS=2*l['N_Sigma']*deposit;RC=l['N_YH']*deposit
        rows.append(dict(name=name,e=e,mA_GeV=massA,Gamma_A_GeV=width,a_mean_decay=decay_a,
                         deposit_R_per_vector=deposit,N_A_band=number,vectors_per_net_charge=number/nx,
                         direct_pair_rest_mass_fraction_if_unprocessed=2*mx*4.8294009076988175e-9*(number/nx)/4.4e-10,
                         R_vector_band_only=R,R_Y=RY,R_Sigma_late_diagnostic=RS,R_cascade_late_diagnostic=RC,
                         combined_late_diagnostic=R+RY+RS+RC,
                         headroom_without_any_scalar_at_R8_2=8.2/deposit/nx,
                         headroom_without_any_scalar_at_R20=20/deposit/nx))
    omega,eq=.8611951308217225,.9233161353726057;Q=200000
    finite=dict(E_per_Q_over_m=eq,chemical_omega_over_m=omega,
                energy_per_charge_excess_fraction=eq/omega-1,
                object_number_chemical_potential_GeV=Q*mx*(eq-omega),
                identity='d[n_L E(Q)] = omega d[n_L Q] + (E-Q omega) d n_L',
                scope='E/Q unequal to dE/dQ is expected for finite objects. A linear lump proxy must include the missing object-number/surface work or use a consistently linear comparator.')
    result=dict(rows=rows,finite_object_thermodynamics=finite,
                primary_comparison='Vector finite band alone, omitting all scalar/cascade and extra release sources.',
                review_thresholds=dict(ke_omega_approx=8.2,ke_EperQ_approx=7.0,ke_with_s_route_bracket=[15.,20.],
                                       source_commit='f7d641a2570f148b67e4ad05bd8acbea5f97c960',
                                       status='Reviewer energy-conserving proxy; not an independently implemented kinetic solver.'))
    write(ROOT/'BUDGET_RESULTS.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
