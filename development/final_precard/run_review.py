"""One-command Issue #9 replay. Private input bytes never enter this repository."""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import platform
import subprocess
import sys
import tempfile
import unittest
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def sha(data): return hashlib.sha256(data).hexdigest()


def json_bytes(value):
    return (json.dumps(value,indent=2,sort_keys=True,default=str,allow_nan=False)+'\n').encode()


def locked(path):
    data=path.read_bytes()
    if sha(data)!=path.with_suffix(path.suffix+'.sha256').read_text().strip():
        raise ValueError('Hash lock mismatch: '+path.name)
    return json.loads(data)


def verified_inputs(path):
    lock=locked(HERE/'INPUT_LOCK.json')
    if sha(path.read_bytes())!=lock['archive_sha256']: raise ValueError('Private input archive changed')
    with zipfile.ZipFile(path) as z:
        names=z.namelist()
        if len(names)!=len(set(names)) or set(names)!=set(lock['members']):
            raise ValueError('Private input archive member inventory changed')
        if any(PurePosixPath(n).is_absolute() or '..' in PurePosixPath(n).parts for n in names):
            raise ValueError('Unsafe archive member path')
        data={n:z.read(n) for n in names}
        if any(sha(b)!=lock['members'][n] for n,b in data.items()):
            raise ValueError('Private input archive member hash mismatch')
    manifest=json.loads(data['INPUTS_MANIFEST.json'])
    if manifest['files']!={n:d for n,d in lock['members'].items() if n!='INPUTS_MANIFEST.json'}:
        raise ValueError('Inner input manifest differs from outer lock')
    for p,blob in lock['canonical_blobs'].items():
        b=data[p]
        if hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()!=blob:
            raise ValueError('Canonical blob identity mismatch: '+p)
    return lock,data


def verify_oracle_freeze():
    frozen=locked(HERE/'ORACLE_FREEZE.json')
    for p,digest in frozen['files'].items():
        if sha((HERE/p).read_bytes())!=digest: raise ValueError('Frozen oracle changed: '+p)
    from . import test_oracles
    suite=unittest.defaultTestLoader.loadTestsFromModule(test_oracles)
    def inventory(s):
        for test in s:
            if isinstance(test,unittest.TestSuite): yield from inventory(test)
            else: yield test.id()
    nodes=sorted(inventory(suite))
    if nodes!=frozen['test_nodes']: raise ValueError('Frozen oracle test inventory drift')
    buffer=io.StringIO(); result=unittest.TextTestRunner(stream=buffer,verbosity=2).run(suite)
    if not result.wasSuccessful(): raise AssertionError(buffer.getvalue())
    return nodes,buffer.getvalue()


def p6_check():
    sys.path.insert(0,str(ROOT/'development'))
    try:
        from integrity import locked_manifest
        from owner_transition import reconciled_manifest
        reconciled_manifest(locked_manifest(ROOT/'development/expected_red_manifest.json'),
                            HERE/'P6_TRANSITION.json',ROOT)
        sys.path.insert(0,str(ROOT/'provenance'))
        import test_doc_sync, expected_red
        suite=unittest.defaultTestLoader.loadTestsFromModule(test_doc_sync)
        def nodes(s):
            for test in s:
                if isinstance(test,unittest.TestSuite): yield from nodes(test)
                else: yield test.id()
        doc_nodes=sorted(nodes(suite))
        buffer=io.StringIO(); result=unittest.TextTestRunner(stream=buffer,verbosity=2).run(suite)
        if not result.wasSuccessful(): raise AssertionError(buffer.getvalue())
        if expected_red._stale_net_cases(): raise ValueError('Stale net prose remains')
        return {'tests':result.testsRun,'pass':result.testsRun,'stale_net_cases':[],
                'test_nodes':doc_nodes,
                'four_doc_scope_verified':True,'P6':'CLOSED_BY_OWNER_RULING_PLUS_NARROW_AMENDMENT'},buffer.getvalue()
    finally:
        sys.path.pop(0); sys.path.pop(0)


def run(inputs,output):
    if platform.python_version()!='3.12.14': raise ValueError('Pinned replay requires Python 3.12.14')
    output.mkdir(parents=True,exist_ok=False)
    lock,payload=verified_inputs(inputs)
    if (HERE/'REVIEW_MANIFEST.json').is_file():
        from .verify_packet import verify
        verify()
    nodes,oracle_log=verify_oracle_freeze()
    p6,p6_log=p6_check()
    (output/'oracle-controls.log').write_text(oracle_log)
    (output/'P6-doc-sync.log').write_text(p6_log)
    # The independent oracle freeze and its exact controls are checked before
    # the first import of ANY supplied/proposed evaluator implementation.
    from development.cr5_review.check import verified_bytes,load,KIT,compare
    from .evaluator_check import prices,baths
    from .nr4_check import regression,weights
    from .weight_boundary import PrevalidatedKit
    with tempfile.TemporaryDirectory(prefix='grut-final-precard-private-') as t:
        t=Path(t)
        for name,data in payload.items():
            path=t/name; path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        cr5_payload,cr5_manifest,_=verified_bytes(t/'cr5-review-bundle.zip')
        for stem in ('nr4_ablation','decoupling_reference'):
            (t/(stem+'.py')).write_bytes(cr5_payload[KIT+stem+'.py'])
        old_l0=load(t/'original/l0_code.py','final_original_l0')
        new_l0=load(t/'proposed/l0_code.py','final_proposed_l0')
        old_hb=load(t/'original/hb_controls.py','final_original_hb')
        new_hb=load(t/'proposed/hb_controls.py','final_proposed_hb')
        prior_l0=load(t/'prior_proposal/l0_code.py','final_prior_proposal_l0')
        price_result=prices(old_l0,new_l0); hb_result=baths(old_hb,new_hb)
        prior_result=prices(old_l0,prior_l0)
        # Keep all newly exposed pre-correction differences, including the
        # lambda alpha-renaming defect and missing explicit prior/table hooks.
        pre_correction=prior_result['mismatches']
        old_binding=sys.modules.get('decoupling_reference')
        try:
            D=load(t/'decoupling_reference.py','decoupling_reference')
            A=load(t/'nr4_ablation.py','final_cr5_nr4')
            pre=load(t/'historical/pre_fix_nr4.py','final_pre_CR4')
            post=load(t/'original/nr4_ablation.py','final_post_CR4')
            nr_result=regression(pre,post,D,A)
            weight_result=weights(A)
            original_rows,original_issues=compare(D,A)
            safe_rows,safe_issues=compare(D,PrevalidatedKit(A))
        finally:
            if old_binding is None: sys.modules.pop('decoupling_reference',None)
            else: sys.modules['decoupling_reference']=old_binding
        if len(original_rows)!=183 or len(safe_rows)!=183:
            raise ValueError('Frozen valid-input comparison inventory changed')
        if original_rows!=safe_rows: raise ValueError('Isolation boundary changed valid-input results')
    mismatches=price_result['mismatches']+hb_result['mismatches']+nr_result['mismatches']+weight_result['mismatches']
    mismatches += [r for r in safe_rows if not r['agrees']]+safe_issues
    status='BLOCKED — semantic/implementation mismatch' if mismatches else 'READY FOR EXTERNAL CHECK'
    record={'status':status,'author_boundary':'SAME_AUTHOR_WORK_EVIDENCE_ONLY',
        'external_checker':None,'external_approval':False,'card_1_authorized':False,
        'cards_consumed':0,'candidate_drafts':0,'bank_approval':False,
        'input_archive_sha256':lock['archive_sha256'],
        'oracle_freeze_sha256':sha((HERE/'ORACLE_FREEZE.json').read_bytes()),
        'oracle_test_nodes':nodes,'oracle_controls_passed':len(nodes),
        'P6':p6,'price_coder':price_result,'HB_controls':hb_result,
        'NR4':nr_result,'CR5_weight_boundary':weight_result,
        'CR5_valid_input_comparisons':safe_rows,'CR5_original_input_contract_findings':original_issues,
        'valid_comparison_count':len(safe_rows),'prior_proposal_mismatches':pre_correction,
        'remaining_semantic_mismatches':mismatches,
        'private_kit_modified':False,'original_direct_calls_cleared':False,
        'review_target':'isolated proposed evaluator plus mandatory PrevalidatedKit boundary',
        'KD1_to_KD5_unbuilt_items':'EXISTING_HOSTILE_DEFAULTS_RETAINED',
        'physical_certificates':'NOT_ESTABLISHED_BY_ABSTRACT_RECEIPTS'}
    inventory=nodes+['P6-doc::'+n for n in p6['test_nodes']]+['price::'+r['fixture'] for r in price_result['rows']]+['HB::'+r['fixture'] for r in hb_result['rows']]
    inventory+=['NR4::'+r['fixture'] for r in nr_result['rows']]+['weights::'+r['fixture'] for r in weight_result['hostile_rows']]
    inventory+=['CR5-183::'+r['fixture'] for r in safe_rows]+['weighted-measure::'+str(i) for i in range(weight_result['additional_valid_measure_comparisons'])]
    if len(inventory)!=len(set(inventory)): raise ValueError('Duplicate review test node')
    (output/'TEST_INVENTORY.json').write_bytes(json_bytes(sorted(inventory)))
    data=json_bytes(record)
    environment={'python':platform.python_version(),'implementation':platform.python_implementation(),
        'platform':platform.platform(),'git_version':subprocess.check_output(['git','--version'],text=True).strip(),
        'actual_checkout_sha':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        'source_references':'all input refs and subject files are hash-pinned; no latest fetch',
        'results_sha256':sha(data),'inventory_sha256':sha((output/'TEST_INVENTORY.json').read_bytes())}
    (output/'RUN_ENVIRONMENT.json').write_bytes(json_bytes(environment))
    # Write a successful result only after the inventory and all provenance
    # stages succeeded. A reporting failure must never leave a READY result.
    (output/'RESULTS.json').write_bytes(data)
    (output/'RESULTS.json.sha256').write_text(sha(data)+'\n')
    print(json.dumps({'status':status,'price_comparisons':price_result['comparison_count'],
        'HB_comparisons':hb_result['comparison_count'],'valid_CR5':len(safe_rows),
        'invalid_entrypoint_checks':weight_result['entrypoint_checks'],
        'additional_valid_measures':weight_result['additional_valid_measure_comparisons'],
        'remaining_mismatches':len(mismatches),'results_sha256':sha(data),
        'test_inventory_count':len(inventory),'external_approval':False,'card_1_authorized':False}))
    return 1 if mismatches else 0


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:
        code=run(args.inputs,args.output)
    except Exception as exc:
        failure={'status':'BLOCKED — '+type(exc).__name__+': '+str(exc),
                 'external_approval':False,'card_1_authorized':False}
        # Do not overwrite a previous run when its output directory already
        # existed. A caller must use a fresh output path for every replay.
        if args.output.is_dir() and not (args.output/'RESULTS.json').exists():
            (args.output/'BLOCKED.json').write_bytes(json_bytes(failure))
        print(json.dumps(failure));code=1
    raise SystemExit(code)
