import json, numpy as np, hashlib, datetime
base='/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl'
refdir='/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019'
mine=json.load(open(f'{base}/reimpl_results.json'))
ser=np.load(f'{base}/series.npz')
ref=json.load(open(f'{refdir}/results/closed_whole_demo.json'))
refm=json.load(open(f'{refdir}/results/marginal_compatibility.json'))
S1=ref['S1_preparation_evolution_measurement']; S2=ref['S2_nonmonotone_relaxation']; S4=ref['S4_initial_correlations']['inf_T_bath_quantum_correlation_Y0X1']; S6=ref['S6_records_and_interventions']; S5=ref['S5_nonunique_local_law']

rows=[]
def cmp(q, m, p, kind='rel', tol=1e-6, note=''):
    m=float(m); p=float(p)
    if kind=='rel':
        d=abs(m-p)/max(abs(p),1e-300); ok=d<=tol; txt=f'rel diff {d:.2e}'
    elif kind=='abs':
        d=abs(m-p); ok=d<=tol; txt=f'abs diff {d:.2e}'
    elif kind=='exact':
        d=abs(m-p); ok=d==0; txt=f'abs diff {d:.2e}'
    rows.append(dict(quantity=q, mine=repr(m), packet=repr(p), agree=bool(ok), note=(txt+('; '+note if note else ''))))
    return ok

# 1
cmp('1. E_tot = Tr(H rho(0))', mine['1_energy']['E_tot_t0'], S1['invariants']['E_tot'], note='packet value printed to 12 significant figures')
cmp('1. max |E_tot(t)-E_tot(0)| on grid', mine['1_energy']['max_abs_drift'], S1['invariants']['max_abs_dE_tot'], kind='abs', tol=1e-12, note='both are double-precision round-off (<1e-13); a 1e-6 RELATIVE criterion is meaningless between two numbers that are themselves numerical zero, so compared as |a-b|<1e-12')
# series check
E_ref=np.array(S1['E_tot']); S_ref=np.array(S1['S_A']); t_ref=np.array(S1['t'])
assert np.allclose(t_ref, ser['grid'])
dE=np.max(np.abs(ser['E_t']-E_ref)); dS=np.max(np.abs(ser['S_A']-S_ref))
cmp('1b. full 401-point E_tot(t) series: max |mine-packet|', dE, 0.0, kind='abs', tol=1e-9, note='series compared pointwise (packet stores 12 sig figs)')
cmp('2b. full 401-point S_A(t) series: max |mine-packet|', dS, 0.0, kind='abs', tol=1e-9, note='series compared pointwise (packet stores 12 sig figs)')
# 2
cmp('2. S_A(t=10)', mine['2_entropy']['S_A_t10'], S_ref[100], note='also packet S5 S_A/10.0/a = %r'%S5['S_A']['10.0']['a'])
cmp('2. S_A(t=40)', mine['2_entropy']['S_A_t40'], S2['S_A_final_t40'])
cmp('2. max_t S_A', mine['2_entropy']['S_A_max'], S2['S_A_max'])
cmp('2. argmax_t S_A', mine['2_entropy']['t_of_max'], S2['t_of_S_A_max'], kind='exact')
# 3
cmp('3. BLP sum of positive increments of D', mine['3_BLP']['sum_positive_increments'], S2['BLP']['backflow_sum_of_positive_increments'])
cmp('3. min_t D(t)', mine['3_BLP']['D_min'], S2['BLP']['D_min'])
cmp('3. argmin_t D(t)', mine['3_BLP']['t_of_min'], S2['BLP']['t_of_D_min'], kind='exact')
cmp('3. D(t=0)', mine['3_BLP']['D_t0'], S2['BLP']['D_t0'])
cmp('3. D(t=40)', float(ser['D_t'][-1]), S2['BLP']['D_t40'])
cmp('3. largest single increase of D', float(np.max(np.diff(ser['D_t']))), S2['BLP']['largest_single_increase'])
# 4
cmp('4. p(Z_0=+1) at t1=3 (before CNOT)', mine['4_record']['p_Z0_plus_t1_before_CNOT'], S6['born_pZ_plus_at_t1'])
cmp('4. p(Z_0=+1) at t2=6, record written', mine['4_record']['p_Z0_plus_t2_with_record'], S6['pZ_plus_A_at_t2_with_record_written'])
cmp('4. p(Z_0=+1) at t2=6, no intervention', mine['4_record']['p_Z0_plus_t2_no_intervention'], S6['pZ_plus_A_at_t2_without_intervention'])
cmp('4. register population |0>_R after write (= p at t1)', mine['4_record']['p_ZR_plus_t2_with_record (record persists)'], S6['register_population_0_after_write'])
# 5
cmp('5. max marginal deviation rho1 vs rho2 (A and Abar)', max(mine['5_initial_correlations']['marginal_A_trace_norm_diff'],mine['5_initial_correlations']['marginal_Abar_trace_norm_diff']), S4['max_marginal_deviation_between_rho1_and_rho2'], kind='exact')
cmp('5. max_t trace distance rho_A(t): rho1 vs rho2', mine['5_initial_correlations']['max_trace_distance_rhoA'], S4['max_trace_distance_rhoA_over_grid'])
cmp('5. argmax_t of that trace distance', mine['5_initial_correlations']['t_of_max'], S4['t_of_max'], kind='exact')
# 6
cmp('6. max Werner p, AB & BC share B (SDP)', mine['6_werner']['p_max_sdp'], refm['C4_werner_thresholds']['chain_AB_BC_max_p'], note='both CLARABEL optima; exact value 2/3 (two-qubit symmetric-extension criterion), solver overshoot ~7e-9')
cmp('6. max Werner p vs analytic 2/3', mine['6_werner']['p_max_sdp'], 2/3, tol=1e-6, note='SDP optimum vs closed form; rel diff 1.07e-8 is solver tolerance')

sha_demo=hashlib.sha256(open(f'{refdir}/code/closed_whole_demo.py','rb').read()).hexdigest()
sha_marg=hashlib.sha256(open(f'{refdir}/code/marginal_compatibility.py','rb').read()).hexdigest()
prov=dict(closed_whole_demo_py_sha256=sha_demo, closed_whole_demo_json_claims=ref['provenance']['script_sha256'], demo_hash_matches=(sha_demo==ref['provenance']['script_sha256']),
          marginal_compatibility_py_sha256=sha_marg, marginal_json_claims=refm['provenance']['script_sha256'], marg_hash_matches=(sha_marg==refm['provenance']['script_sha256']),
          my_reimpl_sha256=hashlib.sha256(open(f'{base}/reimpl.py','rb').read()).hexdigest())
all_agree=all(r['agree'] for r in rows)
json.dump(dict(rows=rows, all_agree=all_agree, provenance=prov), open(f'{base}/comparison.json','w'), indent=2)

# ---- REPORT.md
L=[]
L.append('# Independent re-implementation of the closed-whole demo (from specification only)\n')
L.append(f'Date: {datetime.date.today().isoformat()}  \nPython 3.13 / numpy 2.5.3 / scipy 1.18.1 / cvxpy 1.9.3 (CLARABEL)\n')
L.append('## Code path\n')
L.append(f'- Re-implementation (written before reading the reference script): `{base}/reimpl.py` (sha256 `{prov["my_reimpl_sha256"]}`)')
L.append(f'- Independent expm cross-check of the eigendecomposition evolution: `{base}/xcheck_expm.py`')
L.append(f'- My outputs: `{base}/reimpl_results.json`, time series `{base}/series.npz`')
L.append(f'- Comparison script and output: `{base}/compare.py`, `{base}/comparison.json`')
L.append('- Reference results read AFTER computing mine: `'+refdir+'/results/closed_whole_demo.json`, `'+refdir+'/results/marginal_compatibility.json`')
L.append('- The reference script `closed_whole_demo.py` was NOT read at any point of this task (only its sha256 was taken for provenance).\n')
L.append('Commands:\n```\ncd '+base+'\npython3 reimpl.py          # ~5 s evolution loop + SDP\npython3 xcheck_expm.py     # expm cross-check\npython3 compare.py         # comparison + this report\n```\n')
L.append('## Method\n')
L.append('- Qubit 0 is the leftmost tensor factor; operator on qubit k = I⊗…⊗P_k⊗…⊗I. H built term by term from the spec; asserted Hermitian and equal to 0.5 Z_0 + I⊗H_Abar(7q) + 0.5 X_0X_1.')
L.append('- Bath state exp(−H_Abar)/Z from eigh of the 7-qubit H_Abar (asserted equal to scipy expm version to 1e-12).')
L.append('- Exact evolution U(t)=V e^{−iwt} V† from eigh(H) (256×256); 9-qubit record protocol uses eigh of H⊗I_R (512×512); CNOT = P0_0⊗I + P1_0⊗X_8, asserted CNOT²=I.')
L.append('- Partial trace by reshape+einsum; S_A = −Σλ ln λ; trace norm via SVD; grid t = 0.1·k, k=0..400.')
L.append('- Werner compatibility: single SDP `maximize p s.t. rho_ABC ⪰ 0, Tr=1, Tr_C rho = W(p), Tr_A rho = W(p)`; solution verified (min eigenvalue −1.4e−9, marginal residuals 6e−15) and checked against the closed form p_max = 2/3 from the two-qubit symmetric-extension criterion Tr ρ_B² ≥ Tr ρ_AB² − 4√det ρ_AB (root found by brentq at 0.6666666666666667).\n')
L.append('## My numbers\n```')
L.append(json.dumps(mine, indent=2)); L.append('```\n')
L.append('expm cross-check (independent of eigh): S_A(10)=0.392871183300, S_A(37.5)=0.489338878595, S_A(40)=0.419837916865, D(10)=0.674532015003, D(37.5)=0.543930303492, D(40)=0.674858229496, TD(0.7)=0.10248731291636759, pZ0(3)=0.3549095082640673 — all equal to the eigh numbers to ≤1e-12.\n')
L.append('## Comparison with the packet JSONs\n')
L.append('| quantity | mine | packet | agree | note |\n|---|---|---|---|---|')
for r in rows:
    L.append(f"| {r['quantity']} | {r['mine']} | {r['packet']} | {'yes' if r['agree'] else 'NO'} | {r['note']} |")
L.append('')
L.append(f'**All {len(rows)} comparisons agree: {all_agree}.** The packet stores values to 12 significant figures, so the attainable agreement on scalars is ~1e-12 relative, well inside the 1e-6 criterion. Two quantities are numerical zeros (energy drift 1.8e-14 vs 7.1e-15, marginal deviation 0 vs 0) and are compared absolutely; relative agreement is undefined for them, which is the only "explain" item.\n')
L.append('## Provenance check of the packet\n')
L.append(f'- sha256(code/closed_whole_demo.py) = `{sha_demo}`; JSON claims `{prov["closed_whole_demo_json_claims"]}`; match = {prov["demo_hash_matches"]}')
L.append(f'- sha256(code/marginal_compatibility.py) = `{sha_marg}`; JSON claims `{prov["marginal_json_claims"]}`; match = {prov["marg_hash_matches"]}\n')
L.append('## What the numbers demonstrate (standard QM only, no new physics)\n')
L.append('- Closure forces: E_tot constant to 1e-14 on a mixed global state (unitarity). S_tot/purity constancy and S_A=S_Abar for pure states are the same theorem.')
L.append('- Closure does not force monotone relaxation: S_A rises to 0.489 at t=37.5, sits at 0.420 at t=40, with 13 local maxima; never reaches ln2=0.693.')
L.append('- Closure does not force Markovianity: BLP backflow Σ(ΔD)+ = 1.2035 > 0 (199 of 400 grid steps increase D); a CP-divisible reduced dynamics would give 0.')
L.append('- Records: writing a CNOT record at t=3 changes A\'s later statistics (0.3385 vs 0.4491 at t=6); the register keeps p=0.35491 exactly because [H⊗I_R, Z_R]=0, a property of the supplied H, not of closure.')
L.append('- Initial correlations: rho1=I/256 and rho2=I/256+Y_0X_1/512 have identical A and Abar marginals (deviation exactly 0) but rho_A(t) differ by trace distance 0.1025 at t=0.7; no map on the marginals reproduces both futures (Pechukas).')
L.append('- Marginal compatibility: Werner AB and BC sharing B with identical overlaps admit a common 3-qubit state iff p ≤ 2/3 (SDP 0.66666667, closed form 2/3); above that, overlap agreement is necessary but not sufficient.')
open(f'{base}/REPORT.md','w').write('\n'.join(L)+'\n')
print(json.dumps(dict(rows=rows, all_agree=all_agree, provenance=prov), indent=1))
