"""Review-only weight-contract checks; never modify the supplied Stage-3 kit."""
from fractions import Fraction


def validated_weights(weights):
    """Proposed caller/kit precondition: finite, nonnegative, positive total.

    Conversion and whole-list validation finish before any appearance callback.
    A positive reference measure with zero retained mass is still the existing
    VOID case; a zero reference measure is an invalid normalized-fraction input.
    """
    try:
        values = tuple(Fraction(w) for w in weights)
    except (ValueError, TypeError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('Finite reference weights required') from exc
    if any(w < 0 for w in values) or sum(values) <= 0:
        raise ValueError('Nonnegative reference weights and positive total required')
    return values


def audit(supplied):
    invalid = (
        ('signed-positive-total', [-1, 2]),
        ('negative-second', [2, -1]),
        ('signed-zero-total', [-1, 1]),
        ('zero-reference-measure', [0, 0]),
        ('empty-reference-measure', []),
        ('nan-second', [1, float('nan')]),
        ('positive-infinity-second', [1, float('inf')]),
        ('negative-infinity-second', [1, float('-inf')]),
        ('nonfinite-string-second', [1, 'NaN']),
    )
    rows = []
    for name, weights in invalid:
        try:
            validated_weights(weights)
        except ValueError:
            pass
        else:
            raise AssertionError('Independent proposed precondition accepted invalid input')
        calls = []
        def pipeline(law, i):
            calls.append(i)
            return {'x': i}
        error, result = None, None
        try:
            result = supplied.appearance('fixture', [{'xi': i, 'w': w} for i, w in enumerate(weights)],
                pipeline, lambda image, tol: image['x'] == 1, {'x': 1}, [False] * len(weights))
        except Exception as exc:
            error = type(exc).__name__
        rows.append({'fixture': name, 'weights': [str(w) for w in weights],
                     'required': 'REJECT_BEFORE_APPEARANCE_CALLBACKS',
                     'supplied_exception': error, 'appearance_callbacks': len(calls),
                     'supplied_fraction': None if result is None else str(result['fraction']),
                     'meets_contract': error is not None and not calls})
    # These checks retain positive-total, all-excluded VOID and permit zero atoms.
    assert validated_weights([0, 2]) == (Fraction(0), Fraction(2))
    void = supplied.appearance('fixture', [{'xi': 0, 'w': 1}],
                              lambda *_: (_ for _ in ()).throw(AssertionError('excluded callback')),
                              lambda *_: True, {'x': 1}, [True])
    assert void['fraction'] is None and void['excluded_fraction'] == 1
    return {'owner_contract_url': 'https://github.com/ryangrvr/TestingGRUT/issues/2#issuecomment-6092113130',
            'invalid_input_fixtures': rows, 'violations': sum(not r['meets_contract'] for r in rows),
            'valid_zero_atom_and_positive_measure_VOID_preserved': True,
            'kit_modified': False, 'external_approval': False,
            'status': 'NOT_CLEARED' if any(not r['meets_contract'] for r in rows) else 'CONTRACT_EVIDENCE_ONLY'}
