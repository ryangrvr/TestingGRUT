"""Small exact P1 controls, not a physical-receiver formation solver.

All rooted words, storage configurations and physical charge/field sectors
are retained at the declared small N. Unaddressed fuel qubits are exact
vacuum spectators. Opposite-charge statistics are commuting hard-core slots.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
from math import sqrt
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

JP, JM, KAPPA, EPS = 1/64, 1/96, 1/4, 1/4


@lru_cache(None)
def edges(n):
    return tuple(combinations(range(n), 2))


def fields(n, word, plus, minus):
    """Unique tree Gauss solution in canonical stored label orientation."""
    out = dict.fromkeys(edges(n), 0)
    charge = 0
    for i, j in zip(word, word[1:]):
        charge += int(i in plus)-int(i in minus)
        out[tuple(sorted((i, j)))] = charge if i < j else -charge
    assert charge+int(word[-1] in plus)-int(word[-1] in minus) == 0
    return tuple(out[e] for e in edges(n))


@lru_cache(None)
def matter_states(n, word):
    result = []
    for k in range(len(word)+1):
        for plus in combinations(sorted(word), k):
            for minus in combinations(sorted(word), k):
                ef = fields(n, word, plus, minus)
                if max(map(abs, ef), default=0) <= 1:
                    result.append((plus, minus, ef))
    return tuple(result)


@lru_cache(None)
def occupations(modes, number):
    if modes == 1:
        return ((number,),)
    return tuple((x,)+rest for x in range(number+1)
                 for rest in occupations(modes-1, number-x))


def move(slots, i, j):
    return tuple(sorted(j if x == i else x for x in slots))


def sparse(entries, dim):
    if not entries:
        return coo_matrix((dim, dim), dtype=complex).tocsr()
    rows, cols, vals = zip(*entries)
    return coo_matrix((np.array(vals, complex), (rows, cols)), shape=(dim, dim)).tocsr()


class Model:
    """basis=(word,qmask,receiver tuple,plus positions,minus positions,fuelmask)."""
    def __init__(self, n, matching=((0, 1),), modes=2, fixed_word=None,
                 delta=None, fuel_delta=None, eta=None, lam=None, g=1):
        self.n, self.matching, self.modes, self.g = n, tuple(matching), modes, g
        k = len(matching)
        self.delta = np.array(delta if delta is not None else np.zeros(n))
        self.fdelta = np.array(fuel_delta if fuel_delta is not None else np.zeros(k))
        self.eta = np.array(eta if eta is not None else np.ones(k))
        self.lam = np.array(lam if lam is not None else np.ones(k))
        self.omega = 3+self.fdelta
        if fixed_word is None:
            self.words = [(0,)+p for m in range(1, n+1)
                          for p in permutations(range(1, n), m-1)]
        else:
            self.words = [tuple(fixed_word)]
        self.states, self.ef = [], []
        for w in self.words:
            qr = [sum(1 << i for i in w)] if fixed_word is not None else range(1, 1 << n, 2)
            for q in qr:
                occs = [(0,)*modes] if fixed_word is not None else occupations(modes, q.bit_count()-1)
                for rec, (p, m, ef), f in product(occs, matter_states(n, w), range(1 << k)):
                    self.states.append((w, q, rec, p, m, f)); self.ef.append(ef)
        self.index = {s: i for i, s in enumerate(self.states)}
        self.ef = np.array(self.ef, int)
        self.dim = len(self.states)
        base, fuel, field, particles, density = [], [], [], [], []
        assisted, bare, hop, convert, load = [], [], [], [], []
        self.tags = []

        def put(collection, target, col, val, tag):
            row = self.index[target]
            collection.append((row, col, val)); self.tags.append((row, col, tag))

        for col, (w, q, rec, p, m, f) in enumerate(self.states):
            mask = sum(1 << i for i in w)
            ee = self.ef[col]
            npart = len(p)+len(m)
            fuel_e = sum(om for x, om in enumerate(self.omega) if f >> x & 1)
            qdiag = 0 if fixed_word is not None else (
                16*(q.bit_count()-len(w))+64*(q ^ mask).bit_count()
                +sum(self.delta[i] for i in range(1, n) if q >> i & 1))
            field.append(int(ee @ ee)); particles.append(npart); fuel.append(fuel_e)
            base.append(qdiag+int(ee @ ee)+npart+fuel_e)
            density.append(KAPPA*sum(int(i in p and j in m)+int(i in m and j in p)
                                     for i, j in zip(w, w[1:])))
            if fixed_word is None:
                for x in range(modes-1):
                    for a, b in [(x, x+1), (x+1, x)]:
                        if rec[a]:
                            rr = list(rec); rr[a] -= 1; rr[b] += 1
                            put(assisted, (w, q, tuple(rr), p, m, f), col,
                                -sqrt(rec[a]*(rec[b]+1)), ('receiver',))
                # Directed appends and their adjoints. Fields must remain
                # unchanged; the child is empty. Pop is built only as that
                # exact adjoint, including zero-flux and empty-head gates.
                for v in range(1, n):
                    if v in w or v in p or v in m:
                        continue
                    ww = w+(v,)
                    assert tuple(fields(n, ww, p, m)) == tuple(ee)
                    if not (q >> v & 1):
                        rr = (rec[0]+1,)+rec[1:]
                        target = (ww, q | (1 << v), rr, p, m, f)
                        val = -g/n*sqrt(rec[0]+1)
                        put(assisted, target, col, val, ('append',))
                        assisted.append((col, self.index[target], val))
                    h = w[-1]
                    val = -EPS/sqrt(n)*(int(h in p)+int(h in m))
                    if val:
                        target = (ww, q, rec, p, m, f)
                        put(bare, target, col, val, ('bare',))
                        bare.append((col, self.index[target], val))

            for species, slots, jhop in [(0, p, JP), (1, m, JM)]:
                for old in slots:
                    x = w.index(old)
                    for y in (x-1, x+1):
                        if y < 0 or y >= len(w) or w[y] in slots:
                            continue
                        new = move(slots, old, w[y])
                        pp, mm = (new, m) if species == 0 else (p, new)
                        e2 = fields(n, w, pp, mm)
                        if max(map(abs, e2), default=0) > 1:
                            continue
                        # A hop changes precisely its own physical link.
                        changed = [u for u, (a, b) in enumerate(zip(ee, e2)) if a != b]
                        assert changed == [edges(n).index(tuple(sorted((old, w[y]))))]
                        sign = -1 if old < w[y] else 1
                        assert e2[changed[0]]-ee[changed[0]] == sign*(1 if species == 0 else -1)
                        put(hop, (w, q, rec, pp, mm, f), col, -jhop, ('hop', species))

            bonds = set(zip(w, w[1:]))
            for x, ((i, j), amplitude) in enumerate(zip(matching, self.lam)):
                if (i, j) not in bonds or not (f >> x & 1) or i in p or j in m:
                    continue
                pp, mm = tuple(sorted(p+(i,))), tuple(sorted(m+(j,)))
                e2 = fields(n, w, pp, mm)
                if max(map(abs, e2), default=0) > 1:
                    continue
                link = edges(n).index(tuple(sorted((i, j))))
                assert e2[link]-ee[link] == (1 if i < j else -1)
                target = (w, q, rec, pp, mm, f ^ (1 << x))
                put(convert, target, col, amplitude, ('convert',))
                convert.append((col, self.index[target], amplitude))
            for x, amplitude in enumerate(self.eta):
                target = (w, q, rec, p, m, f ^ (1 << x))
                put(load, target, col, amplitude, ('load',))

        self.field = np.array(field, float)
        self.particles = np.array(particles, float)
        self.fuel = np.array(fuel, float)
        self.fcount = np.array([s[-1].bit_count() for s in self.states], float)
        self.number = self.particles/2
        self.D = diags(base, format='csr', dtype=complex)
        self.G = sparse(assisted, self.dim)
        self.B = sparse(bare, self.dim)
        self.T = sparse(hop, self.dim)
        self.C = sparse(convert, self.dim)
        self.L = sparse(load, self.dim)
        self.V = self.T+diags(density, format='csr')
        self.H0 = self.D+self.G+self.B
        self.HC = self.H0+self.C
        self.Hop = self.H0+self.V
        self.HLrot = self.H0-diags(3*self.fcount)+self.L
        for name in ['H0', 'HC', 'Hop', 'HLrot']:
            h = getattr(self, name)
            assert np.max(np.abs((h-h.getH()).data), initial=0) < 1e-13

    def initial(self, word=None, coherent_receiver=True):
        w = tuple(word if word is not None else range(self.n))
        q = sum(1 << i for i in w)
        recs = [(0,)*self.modes] if len(self.words) == 1 else occupations(self.modes, len(w)-1)
        v = np.zeros(self.dim, complex)
        for z, r in enumerate(recs):
            v[self.index[(w, q, r, (), (), 0)]] = np.exp(.37j*z)/(z+1) if coherent_receiver else int(z == 0)
        return v/np.linalg.norm(v)

    def load(self, state, time, reference=None):
        h = self.HLrot if reference is None else self.reference(self.HLrot, reference)
        return np.exp(-3j*time*self.fcount)*evolve(h, state, time)

    def retained(self, prefix):
        pset = set(prefix)
        return np.array([w[:len(prefix)] == tuple(prefix) and
                         set(p).issubset(pset) and set(m).issubset(pset)
                         for w, q, rec, p, m, f in self.states])

    def reference(self, h, prefix):
        """Common-space compression, with exact reflected-tail modifications."""
        inside = self.retained(prefix)
        co = h.tocoo()
        keep = inside[co.row] == inside[co.col]
        ref = sparse(list(zip(co.row[keep], co.col[keep], co.data[keep])), self.dim)
        # Every bare term is removed only within K. Cut hopping and
        # prefix-breaking assisted transitions already cross K and cancel.
        bb = self.B.tocoo(); sel = inside[bb.row] & inside[bb.col]
        ref -= sparse(list(zip(bb.row[sel], bb.col[sel], bb.data[sel])), self.dim)
        return ref.tocsr()

    def source(self, species, i, j):
        """Gauge-invariant half-bond hopping source (norm <=1/2)."""
        entries = []
        for col, (w, q, r, p, m, f) in enumerate(self.states):
            if (i, j) not in set(zip(w, w[1:])) and (j, i) not in set(zip(w, w[1:])):
                continue
            slots = p if species == 0 else m
            for old, new in [(i, j), (j, i)]:
                if old not in slots or new in slots:
                    continue
                pp, mm = (move(p, old, new), m) if species == 0 else (p, move(m, old, new))
                ef = fields(self.n, w, pp, mm)
                if max(map(abs, ef), default=0) <= 1:
                    entries.append((self.index[(w, q, r, pp, mm, f)], col, .5))
        return sparse(entries, self.dim)

    def number_at(self, species, site):
        return np.array([float(site in s[3+species]) for s in self.states])


def evolve(h, v, time):
    return expm_multiply((-1j*time)*h, v)


def chi(h, psi, source, observable, t):
    u = evolve(h, np.column_stack((psi, source@psi)), t)
    return float(-2*np.imag(np.vdot(u[:, 1], observable*u[:, 0])))


def local_state(delta, eta, lam, tload=np.pi/2, tpump=np.pi/2):
    """Order |0>, |F>, |M>; full pure state, never a Bernoulli substitute."""
    from scipy.linalg import expm
    hr = np.array([[0, eta, 0], [eta, delta, 0], [0, 0, 3]], complex)
    # M is not reachable in loading. Giving it its rotating value is enough
    # for this isolated initial condition; full Model uses exact operator F.
    v = np.exp(-3j*tload*np.array([0, 1, 0]))*(expm(-1j*tload*hr)@np.array([1, 0, 0]))
    hc = np.array([[0, 0, 0], [0, 3+delta, lam], [0, lam, 3]], complex)
    return expm(-1j*tpump*hc)@v


def product_prepared(model, receiver_initial, tload=np.pi/2, tpump=np.pi/2):
    """Exact product pulse reference on a fixed word (or as its embedding)."""
    out = np.zeros(model.dim, complex)
    vectors = [local_state(d, a, b, tload, tpump)
               for d, a, b in zip(model.fdelta, model.eta, model.lam)]
    for col in np.flatnonzero(receiver_initial):
        w, q, r, pp, mm, ff = model.states[col]
        assert not pp and not mm and ff == 0
        for choices in product(range(3), repeat=len(model.matching)):
            p, m, f, amp = [], [], 0, receiver_initial[col]
            for x, ((i, j), z) in enumerate(zip(model.matching, choices)):
                amp *= vectors[x][z]
                if z == 1:
                    f |= 1 << x
                if z == 2:
                    p.append(i); m.append(j)
            out[model.index[(w, q, r, tuple(sorted(p)), tuple(sorted(m)), f)]] += amp
    return out
