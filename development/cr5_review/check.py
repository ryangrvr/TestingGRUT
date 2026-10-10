"""Verify private input hashes first, then compare against a separate oracle.

Private kit files are loaded in a temporary directory, never added to this repo.
Approval registry entries below are synthetic arithmetic fixtures, NOT approved
physical certificates. Separate implementation is not external authorship.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path, PurePosixPath
import sys
import tempfile
import zipfile

from . import oracle as O

BUNDLE_SHA256 = 'f0de023add045314620ee21dddba760de5c3b5e350b4117cb03a46bf2b6eed4b'
MANIFEST_SHA256 = 'ae682955a272c50bd9002c2d833a7cb492bbe082713ae6b9930c789103b40d6b'
RULING_SHA256 = 'd3b5f160735a71ce1017e5e01755f49b548aa9070b534d41c6e54283e6b087a0'
KIT = 'canonical/PROGRAM/STAGE3/kit/'
RULING = 'canonical/PROGRAM/STAGE3/charter_workings/G2_14R_OWNER_RULING.md'


def unique_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('Duplicate manifest key')
            result[key] = value
        return result
    return json.loads(data, object_pairs_hook=pairs)


def verified_bytes(bundle):
    data = Path(bundle).read_bytes()
    if hashlib.sha256(data).hexdigest() != BUNDLE_SHA256:
        raise ValueError('Private review bundle hash mismatch; no code imported')
    with zipfile.ZipFile(bundle) as z:
        names = z.namelist()
        if len(names) != len(set(names)) or any(PurePosixPath(n).is_absolute() or '..' in PurePosixPath(n).parts for n in names):
            raise ValueError('Nonunique or nonlocal bundle member')
        manifest_bytes = z.read('BUNDLE_MANIFEST.json')
        if hashlib.sha256(manifest_bytes).hexdigest() != MANIFEST_SHA256:
            raise ValueError('Bundle manifest does not match frozen review input')
        manifest = unique_json(manifest_bytes)
        if set(names) != set(manifest) | {'BUNDLE_MANIFEST.json'}:
            raise ValueError('Incomplete/extra bundle members')
        payload = {name: z.read(name) for name in manifest}
        for name, digest in manifest.items():
            if hashlib.sha256(payload[name]).hexdigest() != digest:
                raise ValueError('Bundle member hash mismatch: ' + name)
        if hashlib.sha256(payload[RULING]).hexdigest() != RULING_SHA256:
            raise ValueError('Issued normative ruling mismatch')
        cr5 = unique_json(payload['canonical/PROGRAM/RESULTS/WO-004/CR5_MANIFEST.json'])
        for path, digest in cr5['changed_files'].items():
            if manifest.get('canonical/' + path) != digest:
                raise ValueError('CR-5 changed-file manifest disagreement')
        return payload, manifest, cr5


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def compare(D, A):
    rows, issues = [], []
    def check(name, expected, actual):
        rows.append({'fixture': name, 'expected': expected, 'actual': actual, 'agrees': expected == actual})
    def boxes(values):
        return {k: D.Interval(*pair) for k, pair in values.items()}
    binding = D.Binding('full-K', 'context', 'physical-zero', 'exhaustive-domain',
                        'physical-metric', 'finite-chart', 'resolution', 'lock-fiber', ('x',))
    B = O.F
    # Exact grid includes both closed exclusion and strict inclusion boundaries.
    for j, (c, r, e, f, t) in enumerate(itertools.product((-1, 0, 1), (-1, 0, 1), (0, B(1, 4)), (0, B(1, 4)), (0, B(1, 2), 1))):
        a, b, tol = {'x': (c, e)}, {'x': (r, f)}, {'x': t}
        check('interval-' + str(j), O.decide(a, b, tol), D.compare(boxes(a), boxes(b), tol))
    for kind in ('exact', 'oscillation', 'dyadic'):
        probe = {'x': (3, B(1, 16))}
        if kind == 'exact':
            proof = D.ExactZeroProof(binding, 'synthetic-all-branch-proof', boxes(probe))
            expected = O.reference(probe, {'x': 0})
        elif kind == 'oscillation':
            proof = D.OscillationProof(binding, 'synthetic-all-domain-proof', B(1, 4), B(1, 2), boxes(probe), {'x': B(1, 8)})
            expected = O.reference(probe, {'x': B(1, 8)})
        else:
            proof = D.DyadicProof(binding, 'synthetic-annulus-and-tail-proof', 2, boxes(probe), {'x': B(1, 8)}, {'x': B(1, 4)})
            expected = O.reference(probe, {'x': B(3, 8)})
        actual = D.checked_reference('reference', binding, {'reference': proof})
        check('reference-' + kind, expected, {k: (v.center, v.error) for k, v in actual.items()})
    # A receipt is unusable in any changed pipeline/physical context.
    from dataclasses import replace
    registry = {'reference': D.ExactZeroProof(binding, 'synthetic-proof', boxes({'x': (0, 0)})),
                'image': D.ImageProof(binding, 'synthetic-proof', 'point', boxes({'x': (2, 0)}))}
    evidence = D.CornerEvidence('image', 'reference', binding, binding, 'point')
    check('approved-synthetic-corner', 'INCLUDE', D.corner_decision({'x': 2}, evidence, registry, {'x': 1}, full_pipeline='full-K'))
    for field in binding.__dataclass_fields__:
        changed = replace(binding, **{field: ('other',) if field == 'observables' else 'changed'})
        bad = D.CornerEvidence('image', 'reference', changed, binding, 'point')
        check('binding-' + field, 'UNRESOLVED', D.corner_decision({'x': 2}, bad, registry, {'x': 1}, full_pipeline='full-K'))
    for name, ev, reg, sample in (
        ('missing-receipt', replace(evidence, reference_receipt='missing'), registry, {'x': 2}),
        ('wrong-point', replace(evidence, point_digest='elsewhere'), registry, {'x': 2}),
        ('sample-outside-enclosure', evidence, registry, {'x': 3}),
        ('missing-current', evidence, registry, None),
        ('no-approved-registry', evidence, None, {'x': 2}),
    ):
        check(name, 'UNRESOLVED', D.corner_decision(sample, ev, reg, {'x': 1}, full_pipeline='full-K'))
    decisions = ('EXCLUDE', 'INCLUDE', 'UNRESOLVED')
    for ds in itertools.product(decisions, repeat=2):
        for hs in itertools.product((False, True), repeat=2):
            check('nr17-' + repr((ds, hs)), O.nr17(ds, hs), D.nr17_removal_mask(ds, hs))
    for ds in itertools.product(decisions, repeat=2):
        check('q3-' + repr(ds), O.q3(*ds), D.q3_corner_decision(*ds))
    for j, matrix in enumerate((((1, 0), (0, 1)), ((2, 0), (0, B(1, 2))), ((1, 1), (1, -1)))):
        vals = {'x': (2, B(1, 4)), 'y': (-1, B(1, 8))}
        check('transport-' + str(j), O.transport(matrix, vals),
              [(x.center, x.error) for x in D.transport_linear_enclosure(matrix, list(boxes(vals).values()))])
    # Shared full-K mask, missing ablated image retained in denominator, VOID.
    points = [{'xi': i, 'w': 1, 'fiber': 'same'} for i in range(3)]
    registry = {}
    for i in range(3):
        registry['image-' + str(i)] = D.ImageProof(binding, 'synthetic-proof', str(i), boxes({'x': (0 if i == 0 else 2, 0)}))
    registry['reference'] = D.ExactZeroProof(binding, 'synthetic-proof', boxes({'x': (0, 0)}))
    def cert(pt):
        return D.CornerEvidence('image-' + str(pt['xi']), 'reference', binding, binding, str(pt['xi']))
    def pipeline(law, i):
        return {'x': 0 if i == 0 else 2} if law == 'K' else {'x': 2} if law == 'empty' and i == 2 else None
    kwargs = dict(certificates=cert, approved_registry=registry, full_pipeline_digest='full-K')
    result = A.nr4('K', 'empty', 'standard', {'c': 'deleted'}, points, pipeline, None,
                   lambda img, tol: img['x'] == 2, {'x': 1}, 2, **kwargs)
    check('shared-mask', [True, False, False], result['shared_mask'])
    for law, hits in (('K', [False, True, True]), ('K_empty', [False, False, True]), ('K_S', [False]*3), ('K-c', [False]*3)):
        check('shared-measure-' + law, O.measure([1]*3, [True, False, False], hits), result['runs'][law]['fraction'])
    check('undefined-ablation-keeps-denominator', B(1, 2), result['runs']['K_empty']['fraction'])
    check('responsibility', True, result['responsibility']['c'])
    void = A.nr4('K', 'empty', 'standard', {'c': 'deleted'}, points[:1], pipeline, None,
                 lambda img, tol: True, {'x': 1}, 2, **kwargs)
    check('nr4-void-fractions', [None]*4, [r['fraction'] for r in void['runs'].values()])
    check('nr4-void-responsibility', 'VOID', void['responsibility']['c'])
    # Hostile inputs are deliberately NOT silently repaired in the supplied kit.
    hostile = [{'xi': 0, 'w': -1}, {'xi': 1, 'w': 2}]
    try:
        bad = A.appearance('K', hostile, lambda law, i: {'x': i}, lambda img, tol: img['x'] == 1,
                           {'x': 1}, [False, False])
    except (ValueError, D.CertificateError):
        pass
    else:
        issues.append({'kind': 'UNVALIDATED_REFERENCE_WEIGHTS', 'expected': 'reject negative reference-measure weights',
                       'actual_fraction': bad['fraction'], 'severity': 'INPUT_CONTRACT_REVIEW_REQUIRED',
                       'scope': 'invalid input; not evidence that valid certified arithmetic fails'})
    return rows, issues


def run(bundle, output):
    payload, manifest, governance = verified_bytes(bundle)  # FIRST, before any code loading.
    output.mkdir(parents=True, exist_ok=False)
    checkpoint = {'oracle_sha256': hashlib.sha256(Path(O.__file__).read_bytes()).hexdigest(),
                  'controls_passed': O.controls(), 'supplied_code_imported': False,
                  'normative_ruling_sha256': RULING_SHA256, 'verified_payload_files': len(manifest)}
    (output / 'independent-ready.json').write_text(json.dumps(checkpoint, indent=2) + '\n')
    # Compare only after the separate implementation and its controls exist.
    with tempfile.TemporaryDirectory(prefix='grut-private-cr5-') as tmp:
        tmp = Path(tmp)
        for stem in ('decoupling_reference', 'nr4_ablation'):
            (tmp / (stem + '.py')).write_bytes(payload[KIT + stem + '.py'])
        previous = sys.modules.get('decoupling_reference')
        try:
            D = load(tmp / 'decoupling_reference.py', 'decoupling_reference')
            A = load(tmp / 'nr4_ablation.py', 'review_supplied_nr4')
            rows, issues = compare(D, A)
            from .weight_contract import audit
            weight_contract = audit(A)
        finally:
            if previous is None:
                sys.modules.pop('decoupling_reference', None)
            else:
                sys.modules['decoupling_reference'] = previous
            sys.modules.pop('review_supplied_nr4', None)
    discrepancies = [r for r in rows if not r['agrees']]
    result = {'scope': 'SEPARATE_IMPLEMENTATION_CHECK_NOT_EXTERNAL_APPROVAL',
              'bundle_sha256': BUNDLE_SHA256, 'normative_ruling_sha256': RULING_SHA256,
              'hashes_verified_before_import': True, 'payload_members_verified': len(manifest),
              'comparison_count': len(rows), 'discrepancies': discrepancies,
              'input_contract_findings': issues, 'comparisons': rows,
              'owner_weight_contract': weight_contract,
              'valid_input_arithmetic': 'PASS' if not discrepancies else 'FAIL',
              'broader_evaluator': 'NOT_IN_THIS_BUNDLE_NOT_REVIEWED',
              'certificate_theorems_and_physical_locks': 'NOT_CERTIFIED_BY_SOFTWARE_FIXTURES',
              'private_source_committed': False, 'canonical_publication': governance['canonical_publication'],
              'canonical_re_freeze': governance['canonical_re_freeze'],
              'external_approval': False, 'card_1_authorized': False,
              'verdict': 'REVIEW_REQUIRED' if not discrepancies else 'DISCREPANCY'}
    (output / 'comparison.json').write_text(json.dumps(result, indent=2, default=str) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'comparisons'}, default=str), flush=True)
    return 1 if discrepancies else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.bundle, args.output))
