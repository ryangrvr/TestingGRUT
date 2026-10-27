"""Population-independent decay clocks; no recomputed source normalization."""
from timing_common import *
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

def energy_tail(tau,gam,C,offset):
    u=gam*(tau-offset)
    assert u>0
    return float(np.sqrt(2*C/gam)*np.exp(-gam*offset)*gamma(1.5)*gammaincc(1.5,u))

class Clock:
    def __init__(self,times,a,C,offset,gam,order=8):
        self.tm=float(times[-1]);self.C=C;self.offset=offset;self.gam=gam
        self.times=np.asarray(times);self.order=order
        self.af=PchipInterpolator(times,a,extrapolate=False)
        self.early=self.energy_early(self.tm)
        self.total=self.early+energy_tail(self.tm,gam,C,offset)
    def a_of_tau(self,tau):
        return float(self.af(tau)) if tau<=self.tm else float(np.sqrt(2*self.C*(tau-self.offset)))
    def energy_early(self,tau):
        # Integrate each cubic interpolation interval explicitly. This avoids
        # adaptive quadrature warnings at thousands of interpolation knots.
        edges=np.unique(np.r_[0,self.times[(self.times>0)&(self.times<tau)],tau])
        nodes,weights=np.polynomial.legendre.leggauss(self.order)
        half=np.diff(edges)/2;mid=(edges[1:]+edges[:-1])/2
        tt=mid[:,None]+half[:,None]*nodes[None,:]
        value=self.gam*self.af(tt)*np.exp(-self.gam*tt)
        return float(np.sum(half[:,None]*weights[None,:]*value))
    def number_fraction(self,tau):return float(-np.expm1(-self.gam*tau))
    def energy_fraction(self,tau):
        value=self.energy_early(tau) if tau<=self.tm else self.early+energy_tail(self.tm,self.gam,self.C,self.offset)-energy_tail(tau,self.gam,self.C,self.offset)
        return float(value/self.total)
    def quantile(self,fraction,kind='number'):
        if kind=='number':tau=-np.log1p(-fraction)/self.gam
        else:
            u=brentq(lambda u:self.energy_fraction(u/self.gam)-fraction,1e-8,30,xtol=1e-13)
            tau=u/self.gam
        return dict(fraction=fraction,tau=float(tau),a=self.a_of_tau(tau))

def main():
    result=read('RESPONSE_RESULTS.json');out=[]
    widths={'YL':K.SRC['mediator_debris']['YL_three_body']['parents']['Gamma_L_GeV'],
        'A_e1e-7':K.SRC['vector_sources']['e1e-7_band']['parent']['Gamma_A_GeV'],
        'A_e1e-6':K.SRC['vector_sources']['e1e-6_band']['parent']['Gamma_A_GeV']}
    import reheating as R
    assert R.H0==K.H_I
    mL=R.H0*np.sqrt(R.M**2-np.sqrt(2)*R.MU*R.HH)
    matching={'H_i_relative':R.H0/K.H_I-1,
        'YL_width_relative':3*R.KAP**2*mL/(256*np.pi**3)/widths['YL']-1}
    for e,label in [(1e-7,'A_e1e-7'),(1e-6,'A_e1e-6')]:
        mA=6*e*R.HH*R.H0
        matching[label+'_width_relative']=e*e*mA/(48*np.pi)*(1-4*K.MX**2/mA**2)**1.5/widths[label]-1
    assert max(abs(v) for v in matching.values())<1e-12
    for row in result['rows']:
        if not row['tag'].startswith('quiet'):continue
        z=np.load(ROOT/(row['tag']+'.npz'));times=z['times'];a=np.exp(z['bg'][:,15]);C=row['late_C'];offset=row['matched_radiation_tau_offset']
        clocks={}
        for name,width in widths.items():
            clock=Clock(times,a,C,offset,width/K.H_I);crossings=[]
            for event in row['crossings']:
                for kind in ['amplitude_10percent','amplitude_unity','local_charge_unity']:
                    crossing=event[kind]
                    if crossing is None:continue
                    tt=crossing['tau_bracket'][1]
                    crossings.append(dict(phase_seed=event['phase_seed_radians'],kind=kind,tau=tt,a=clock.a_of_tau(tt),
                        decayed_number_fraction_rest_clock=clock.number_fraction(tt),
                        cumulative_conformal_energy_fraction_prescribed_source=clock.energy_fraction(tt)))
            clocks[name]=dict(Gamma_GeV=width,gamma_over_Hi=width/K.H_I,
                number_quantiles=[clock.quantile(f,'number') for f in [.01,.1,.5,.9,.99]],
                energy_quantiles=[clock.quantile(f,'energy') for f in [.1,.5,.9]],
                mean_a_per_decay=clock.total,at_response_crossings=crossings,
                late_RD_asymptotic_mean_a=float(np.sqrt(np.pi*C/2/(width/K.H_I))),
                old_instant_RD_mean_a=float(np.sqrt(np.pi/2/(width/K.H_I))),
                caveat='Pinned late fixed width, parents available at tau0, nonrelativistic clock. Early rest-clock number is an upper bound with relativistic time dilation. No production population, daughter thermalization or actual deposited heat computed.')
        aq=(abs(row['late_NX_normalization_through80'])*K.H_I**3/K.N_Q0)**(1/3)
        kinematic=[]
        for k in [2.13,5.8]:
            kinematic.append(dict(k_over_Hi=k,Q_in_full_wavelength_cube=abs(row['late_NX_normalization_through80'])*(2*np.pi/k)**3,
                p_over_mX_at_reference_density=k*K.H_I/(aq*K.MX),
                scope='Kinematic charge in a cube of side one wavelength; not a fragment distribution, nonlinear object charge or a formation rule.'))
        out.append(dict(tag=row['tag'],late_C=C,matched_radiation_tau_offset=offset,sources=clocks,
            a_at_old_reference_charge_density=aq,reference_density_is_not_formation_prediction=True,
            wavelength_diagnostics=kinematic))
    # Quantify the review's neutral-reservoir note without assigning quanta an origin.
    cw=ROOT.parent/'grut_core_wave_audit_20261026'
    state=json.loads((cw/'FINAL_AUDIT.json').read_text())['cases']['mixed_N96']['whole_box']
    V0=75*48**3;EX0=67.5*48**3
    neutral=dict(final_Ephi_over_initial_V0=state['E_phi']/V0,
        net_charged_energy_gain_over_V0=(state['E_X']-EX0)/V0,
        total_energy_drift_over_V0=(state['E_total']-V0-EX0)/V0,
        scope='Energy comparison under the stated interaction convention, not tagging individual quanta or a universal clump-energy partition.')
    save('SOURCE_CLOCKS.json',dict(histories=out,neutral_reservoir_note=neutral,width_matching_same_couplings=matching,
        parent_widths_reused=True,production_spectra_and_normalizations_recomputed=False,
        physical_f_sym=None,scope='Conditional growth/source chronology. The old R normalization is not transplanted onto these reheating histories.'))
    for row in out:
        print('CLOCKS',row['tag'],'C',row['late_C'],'a_reference',row['a_at_old_reference_charge_density'],
              [(key,v['number_quantiles'][2]['a']) for key,v in row['sources'].items()],flush=True)

if __name__=='__main__':main()
