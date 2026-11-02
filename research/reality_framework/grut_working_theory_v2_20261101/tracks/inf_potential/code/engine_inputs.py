#!/usr/bin/env python3
"""inf_potential (Claude, own code): engine inputs for rpinf.c MODE 4 at the S1 slice (a = 1, H = H_i), per benchmark and H_i.
G = (10i - mu)/phihat_I(S1) so that mu_eff = 10i at S1; for the hilltop the same holder written in the field measured from the
top of the hill has the complex constant mu_b = mu - G phihat_min.  Output ../results/engine_inputs.json."""
import json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inf_common as ic
RES = Path(__file__).resolve().parent.parent / 'results'
mu = 7.159456159
out = {}
for name in ('HT6t', 'T1', 'E1'):
    d = json.load(open(RES / ('bench_%s.json' % name)))
    for r in d['by_Hi']:
        Hi = r['H_i_GeV']; e2 = (Hi / ic.MPL_GEV) ** 2
        m = ic.Model({'T1': 'T', 'E1': 'E', 'HT6t': 'HT'}[name], **r['params'])
        I0 = r['phiI_S1_MP']; G = (10j - mu) / (I0 / math.sqrt(e2))
        rec = dict(eta2=e2, b_cap=r['m_I_over_Hi'], RPINF_KIND={'T1': 1, 'E1': 2, 'HT6t': 3}[name], RPINF_AMP=r['V0_MP4'] / e2,
                   RPINF_P1=m.vev if name == 'HT6t' else m.M, RPINF_PHIMIN=m.phi_min(), RPINF_I0=I0,
                   RPINF_SIGN=1 if r['phidot_S1_over_HiMP'] > 0 else -1, Ip0=r['phidot_S1_over_HiMP'], G=[G.real, G.imag], g_I=abs(G),
                   alpha_I=math.atan2(G.imag, G.real))
        if name == 'HT6t':
            rec.update(RPINF_P2=m.pw, RPINF_P3=m.beta, RPINF_FTOP=m.ftop)
            mub = mu - G * (m.phi_min() / math.sqrt(e2)); rec['mu_bare_in_phi_frame'] = [mub.real, mub.imag]
        out['%s_%.3g' % (name, Hi)] = rec
json.dump(dict(note=__doc__, inputs=out), open(RES / 'engine_inputs.json', 'w'), indent=1)
print(json.dumps(out, indent=0)[:400])
