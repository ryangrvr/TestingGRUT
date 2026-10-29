"""Count each supported core reaction along the unchanged kinetic history."""
from support import *
import transport as tr
from scipy.integrate import solve_ivp
import argparse, time

CHANNELS = ['X_to_Xbar', 'Xbar_to_X', 'X_to_h_then_pair', 'Xbar_to_h_then_pair']
# Columns follow CHANNELS; rows are Nplus, Nminus, Qcore.
STOICH = np.array([[-1, 1, 0, 1], [1, -1, 1, 0], [2, -2, 1, -1.]])

def run(bins, rtol, tag):
    started=time.monotonic()
    e=tr.grid(0,8,'translated',bins)
    model=tr.Model(e,'full',12); n=len(e); dim=2*n
    y0,omitted=tr.initial(e)
    rates=model.core_rates
    rows=np.vstack((np.r_[np.ones(n),np.zeros(n)],np.r_[np.zeros(n),np.ones(n)]))
    operator_stoich_error=float(np.max(abs(rows@model.L-STOICH[:2]@rates)))
    def rhs(t,y):
        dy,lost=model.rates(y[:dim])
        return np.r_[dy,lost,rates@y[:dim]]
    def jac(t,y):
        j=np.zeros((dim+5,dim+5));j[:dim,:dim]=model.rates(y[:dim],True)
        j[dim+1:,:dim]=rates
        return j
    result=solve_ivp(rhs,(0,1000),np.r_[y0,np.zeros(5)],method='BDF',jac=jac,
                     rtol=rtol,atol=1e-9,dense_output=True)
    assert result.success,result.message
    reference=np.load(KIN/f'translated_full_{bins}.npz')
    times=reference['tau']; y=result.sol(times); counts=y[dim+1:]
    totals=rows@y[:dim]; changes=totals-totals[:,0,None]
    ledger_error=float(np.max(abs(changes-STOICH[:2]@counts)))
    pop_difference=float(np.max(abs(y[:dim]-reference['Y'][:dim])))
    endpoint_relative=np.max(abs(totals[:,-1]/(rows@reference['Y'][:dim,-1])-1))
    q=fixed.QREF-totals[0]+totals[1]
    invariant=np.r_[e-fixed.OMEGA,e+fixed.OMEGA]@y[:dim]
    counts_final=counts[:,-1]; reverse,forward,positive_h,negative_h=counts_final
    net=forward-reverse
    decomposition={'net_forward_exchange':net,
       'core_charge_change_exchange':-2*net,
       'core_charge_change_neutral':positive_h-negative_h,
       'total_neutral_events':positive_h+negative_h,
       'core_charge_change_total':q[-1]-q[0],
       'net_antiparticle_loss':totals[1,0]-totals[1,-1],
       'gas_total_number_increase':float(changes[:,-1].sum())}
    out={'configuration':{'shift_bins':bins,'angles':12,'tau_final':1000,'rtol':rtol,'atol':1e-9},
      'channel_order':CHANNELS,'stoichiometry_rows_Nplus_Nminus_Qcore':STOICH,
      'integrated_supported_events':dict(zip(CHANNELS,counts_final)),
      'decomposition':decomposition,'initial':{'Nplus':totals[0,0],'Nminus':totals[1,0],'Qcore':q[0]},
      'final':{'Nplus':totals[0,-1],'Nminus':totals[1,-1],'Qcore':q[-1]},
      'checks':{'operator_stoichiometry_error':operator_stoich_error,'integrated_stoichiometry_error':ledger_error,
        'maximum_bin_population_difference_from_parent':pop_difference,
        'endpoint_species_relative_difference_from_parent':float(endpoint_relative),
        'invariant_relative_drift':float(np.max(abs(invariant/invariant[0]-1))),
        'minimum_population':float(y[:dim].min()),'omitted_events':float(y[dim,-1]),
        'initial_omitted_support':omitted},
      'solver':{'success':bool(result.success),'nfev':result.nfev,'steps':len(result.t),'seconds':time.monotonic()-started}}
    np.savez_compressed(ROOT/(tag+'.npz'),tau=times,counts=counts,N=totals,Qcore=q,
                        gas_invariant=invariant,omitted_events=y[dim])
    save(tag+'.json',out);print(json.dumps(out,default=lambda x:x.tolist() if hasattr(x,'tolist') else x),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--bins',type=int,default=96)
    ap.add_argument('--rtol',type=float,default=2e-7);ap.add_argument('--tag',default='channels96')
    a=ap.parse_args();run(a.bins,a.rtol,a.tag)
