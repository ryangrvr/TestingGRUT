"""Unbiased-by-phi check: select high current first, then inspect its environment."""
from pilot import *

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];rows=[]
    for tag,r in runs.items():
        N,L,t,a=snapshot(r['snapshots'][-1]);j=np.asarray(a[0])*a[4]-np.asarray(a[1])*a[3];mass2=np.asarray(a[2])**2/C.NU**2
        high=j>=np.quantile(j,.9)
        corr=float(np.mean((j-j.mean())*(mass2-mass2.mean()))/(j.std()*mass2.std()))
        rows.append(dict(tag=tag,pearson_current_vs_phi_squared=corr,
            mean_mass_squared_over_m2=float(mass2.mean()),high_current_decile_mean_mass_squared_over_m2=float(mass2[high].mean()),
            depleted_fraction_all=float(np.mean(mass2<.25)),depleted_fraction_high_current_decile=float(np.mean(mass2[high]<.25))))
    save('ASSOCIATION.json',dict(rows=rows,
        selection='Top10% current cells, selected without looking at phi. phi²/v² is the binding-field contribution g²phi²/m_X²; it excludes the quartic and other medium contributions to fluctuation masses.',
        scope='Spatial association in these supplied cold histories, not an independent causal intervention, a particle population or a universal law.'))
    print([(r['tag'],r['pearson_current_vs_phi_squared'],r['high_current_decile_mean_mass_squared_over_m2']) for r in rows],flush=True)

if __name__=='__main__':main()
