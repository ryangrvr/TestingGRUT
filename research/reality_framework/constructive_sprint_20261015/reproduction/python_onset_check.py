"""Fresh cross-language signed-path check for one published order-22 target.

This evaluates the same monotone-path recurrence in Python integers and uses
the frozen Python fermion-hop implementation. It is a language/implementation
control, not a separately invented physical model or a second proof.
"""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import json
import sys
import time

P = Path(__file__).resolve().parent
ROOT = P.parents[1]
sys.path.insert(0, str(ROOT / "cp12_ground_connectivity"))
from exact_destination import gauss, hop, sea

def main():
    start = time.monotonic()
    n = 5
    es = list(combinations(range(n), 2))
    holes = {(0, 1), (1, 2)}
    graph = (1 << len(es)) - 1
    fields = []
    charged = set()
    for i, j in es:
        tail, head = (i, j) if (j-i) % n in (1, 2) else (j, i)
        is_hole = (tail, head) in holes
        if is_hole:
            charged.add(tail)
        fields.extend([1 + int(is_hole) if tail == i else -1,
                       1 + int(is_hole) if tail == j else -1])
    bits = sum(1 << (n+e) for e, pair in enumerate(es) if pair not in holes) + sum(1 << i for i in charged)
    final = (graph, bits, tuple(fields))
    gauss(n, final)
    double = sum(1 << j for j, x in enumerate(fields) if abs(x) == 2)
    v = {(0, 0, sea(n, graph)[1]): 1}
    counts = [1]
    order = sum(map(abs, fields))
    for depth in range(order):
        nxt = defaultdict(int)
        for (once, twice, fermions), amplitude in v.items():
            for j, field in enumerate(fields):
                bit = 1 << j
                if once & bit:
                    if not double & bit or twice & bit:
                        continue
                    new_once, new_twice = once, twice | bit
                else:
                    new_once, new_twice = once | bit, twice
                edge, leg = divmod(j, 2)
                vertex = es[edge][leg]
                dst, src = (vertex, n+edge) if field > 0 else (n+edge, vertex)
                moved = hop(fermions, dst, src)
                if moved is not None:
                    new_bits, sign = moved
                    nxt[new_once, new_twice, new_bits] += amplitude * sign * (1 if leg == 0 else -1)
        v = {key: value for key, value in nxt.items() if value}
        counts.append(len(v))
    amplitude = v.get(((1 << len(fields))-1, double, bits), 0)
    out = dict(status="PASS" if amplitude == -49583226880 else "FAIL",
        exact_signed_amplitude=str(amplitude), expected_signed_amplitude="-49583226880",
        minimum_hops=order, dynamic_program_layer_sizes=counts,
        target=dict(graph=graph, fermion_bits=bits, electric_fields=fields),
        seconds=time.monotonic()-start,
        scope="Fresh Python cross-language control of one exact shortest-path amplitude; no full time evolution.")
    (P / "PYTHON_ONSET_CHECK.json").write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps(out, indent=2))
    raise SystemExit(out["status"] != "PASS")

if __name__ == "__main__":
    main()
