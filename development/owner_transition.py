"""Apply the exact Issue #2 owner reconciliation to an observational comparison.

The original manifest stays locked and unchanged. This is not a general refresh
or permission to change outcomes, bank decisions, seals, or scientific inputs.
"""
from copy import deepcopy
import hashlib
from pathlib import Path
import subprocess

from integrity import json_read, source_identity

OWNER_URL = 'https://github.com/ryangrvr/TestingGRUT/issues/2#issuecomment-6092113130'
POINTER = ('test_prereg_immutable.py::TestBlindSafe::'
           'test_no_sealed_prereg_points_outward_at_its_own_context')
ADDED = {f'PREREG_TERMINATION_V{v}_2026-08-10.txt -> RESULT_TERMINATION_events.txt :: {kind}':
         'P2-TERMINATION-EVENTLOG' for v in (3, 4) for kind in (
             'a sigma value', 'a sigma value stated as a quantity rather than as a bound')}
P1A = 'P1A-EDGE-REPRESENTATION'
PATHS = {'provenance/expected_red.py', 'provenance/OPEN_PASSES.txt'}


def reconciled_manifest(base, path, root):
    path, root = Path(path), Path(root)
    payload = path.read_bytes()
    lock = path.with_suffix(path.suffix + '.sha256').read_text().strip()
    if hashlib.sha256(payload).hexdigest() != lock:
        raise ValueError('Owner transition hash lock mismatch')
    transition = json_read(path)
    if transition['owner_direction_url'] != OWNER_URL or transition['scope'] != 'WORKBENCH_ONLY_NO_SCIENTIFIC_ACCEPTANCE':
        raise ValueError('Unexpected owner reconciliation authority/scope')
    original_path = root / 'development/expected_red_manifest.json'
    if transition['original_manifest_sha256'] != hashlib.sha256(original_path.read_bytes()).hexdigest():
        raise ValueError('Owner transition references a different original manifest')
    if transition['base_protected_inputs'] != base['protected_inputs']:
        raise ValueError('Owner transition changes its base protected inputs')
    direction = root / 'development/owner_reconciliation/OWNER_DIRECTION.md'
    if hashlib.sha256(direction.read_bytes()).hexdigest() != transition['owner_direction_file_sha256']:
        raise ValueError('Recorded owner direction changed')
    if set(transition['changed_protected_paths']) != PATHS:
        raise ValueError('Owner transition touches an unauthorized protected path')
    for p, digest in transition['changed_protected_paths'].items():
        if hashlib.sha256((root / p).read_bytes()).hexdigest() != digest:
            raise ValueError('Owner reconciliation source drift: ' + p)
    changed = subprocess.check_output(['git', 'diff', base['source_commit'], '--name-only'],
                                      cwd=root, text=True).splitlines()
    protected = {p for p in changed if p != 'AGENT_COORDINATION.md'
                 and not p.startswith(('development/', '.github/workflows/'))}
    identity = source_identity(root)
    if protected != PATHS or identity['untracked_protected_paths']:
        raise ValueError('Additional protected changes outside owner reconciliation')
    observed = identity['protected_inputs']
    if (observed != transition['reconciled_protected_inputs']
            or observed['tracked_file_count'] != base['protected_inputs']['tracked_file_count']):
        raise ValueError('Reconciled protected-input inventory mismatch')
    expected = deepcopy(base)
    entry = expected['declarations'][POINTER]
    if any(k in entry['cases'] for k in ADDED):
        raise ValueError('P2 cases already present in original manifest')
    entry['cases'].update(ADDED)
    # Live enumerators already saw these cases in the original snapshot. Only
    # the declared ownership changes; neither scanner nor live case set changes.
    if not set(ADDED) <= set(entry['live_cases']):
        raise ValueError('P2 scope would hide a newly changed live enumerator')
    expected['open_passes'][P1A]['status'] = 'CLOSED'
    expected['protected_inputs'] = transition['reconciled_protected_inputs']
    if transition['declarations'] != expected['declarations'] or transition['open_passes'] != expected['open_passes']:
        raise ValueError('Unapproved declaration or pass disposition in transition')
    # No table from the transition may overwrite these independent observations.
    if set(transition) != {'schema_version', 'scope', 'owner_direction_url',
                          'owner_direction_file_sha256', 'original_manifest_sha256',
                          'base_protected_inputs', 'reconciled_protected_inputs',
                          'changed_protected_paths', 'declarations', 'open_passes'}:
        raise ValueError('Unexpected transition fields; scientific/outcome refresh forbidden')
    return expected
