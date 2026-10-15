"""Unregulated signed-amplitude onset certificates for supplied charged cycles."""
from fractions import Fraction
from math import factorial
import json

from run import HERE, initial
from exact_destination import powers
from frozen_model import graph_basis


def main():
    g = initial(4, 1)[0][0]
    # Any zero edge fixes the sole cycle flow.  Gauss constraints then bound each
    # remaining half-leg by four, making this a complete finite target catalog.
    target = [s for s in graph_basis(4, g, 4)
              if any(s[2][2*e] == s[2][2*e+1] == 0 for e in range(6) if g >> e & 1)]
    maximum_target_field = max(abs(x) for s in target for x in s[2])
    assert maximum_target_field == 2
    distance_rows = []
    for L in range(1, 13):
        s, _ = initial(4, L)
        d = min(sum(abs(a-b) for a,b in zip(s[2], t[2])) for t in target)
        assert d == 8*L-5
        distance_rows.append(dict(flux=L, minimum_electric_L1_distance_to_zero_edge=d))
    # For L>=2 every target half-leg lies between the initial positive and
    # negative magnitudes. Each of the eight nonzero initial legs therefore
    # contributes slope +1, so the finite check at L=2 proves the whole tail.
    rows = []
    for L in [1,2,3,4]:
        s, sign = initial(4, L)
        order = 8*L-4
        scale, vectors = powers(4, {s: sign}, order)
        first_zero = None
        first_change = None
        for k, v in enumerate(vectors):
            zero = {t:a for t,a in v.items()
                    if any(t[2][2*e] == t[2][2*e+1] == 0 for e in range(6) if g >> e & 1)}
            changed = {t:a for t,a in v.items() if t[0] != g}
            if zero and first_zero is None:
                first_zero = k
            if changed and first_change is None:
                first_change = k
        assert first_zero == 8*L-5
        assert first_change == 8*L-4
        changed = {t:a for t,a in vectors[first_change].items() if t[0] != g}
        coefficient = Fraction(sum(a*a for a in changed.values()),
                               scale**(2*first_change)*factorial(first_change)**2)
        rows.append(dict(flux=L, first_zero_edge_H_power=first_zero,
            first_changed_graph_H_power=first_change,
            graph_change_probability_leading_power=2*first_change,
            signed_leading_target_configurations=len(changed),
            graph_change_probability_coefficient=str(coefficient),
            coefficient_float=float(coefficient),
            reachable_support_counts=[len(v) for v in vectors]))
    out = dict(status="PASS", zero_edge_target_configurations=len(target),
        target_catalog_completeness="A zero half-leg pair fixes the only cycle-flow variable. Following the alternating vertex and edge Gauss equations around four vertices changes the field by at most one at each step, so |E|<=4 contains every target; finite Fock occupancy is exhaustively enumerated.",
        maximum_absolute_target_field=maximum_target_field,
        all_integer_flux_distance_formula="For every integer L>=1: d_min=8L-5. L=1 checked directly. All target fields satisfy |E|<=2; for L>=2 each target's L1 distance is an affine function of slope8, so L=2 fixes the entire remaining range.",
        minimum_field_distances=distance_rows, exact_unregulated_onsets=rows,
        interpretation="Prepared charged four-cycles only. Powers use the complete frozen Hamiltonian, integer signed amplitudes, no electric cutoff and no path-count substitution. These calculations certify nonzero leading loss coefficients, not macroscopic survival or autonomous preparation.")
    (HERE / "EXACT_ONSET.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
