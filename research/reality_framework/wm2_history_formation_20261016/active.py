"""History-exact bosonic mode compression with an accumulated error account.

Only the receiver one-particle subspace is compressed. Every type word,
reverse transition and charge-allowed occupation in that subspace is kept.
The omitted state norm is counted, not declared physically inaccessible.
"""
from functools import lru_cache
from itertools import combinations_with_replacement
from math import ceil, comb, exp, factorial, sqrt
import time

import numpy as np


def type_words(counts):
    result = []

    def visit(w):
        result.append(w)
        for a, n in enumerate(counts):
            if w.count(a) < n:
                visit(w + (a,))

    visit(())
    return sorted(result, key=lambda w: (len(w), w))


class Layout:
    def __init__(self, counts, values, rank, state_cap=2_000_000):
        self.counts, self.values, self.rank = tuple(counts), tuple(values), rank
        self.M = sum(counts)
        self.words = type_words(counts)
        estimated = sum((comb(rank+len(w)-1, len(w)) if rank else int(not w)) for w in self.words)
        if estimated > state_cap:
            raise ResourceLimit(f"coefficient dimension {estimated} exceeds {state_cap}")
        self.word_index = {w: j for j, w in enumerate(self.words)}
        self.occ = [list(combinations_with_replacement(range(rank), k))
                    for k in range(self.M + 1)]
        self.index = [{o: j for j, o in enumerate(a)} for a in self.occ]
        self.nwords = [sum(len(w) == k for w in self.words) for k in range(self.M + 1)]
        self.koff = np.r_[0, np.cumsum([n * len(o) for n, o in zip(self.nwords, self.occ)])]
        self.size = int(self.koff[-1])
        if self.size > state_cap:
            raise ResourceLimit(f"coefficient dimension {self.size} exceeds {state_cap}")
        self.offset, self.energy = [], np.empty(self.size)
        offset = 0
        for w in self.words:
            self.offset.append(offset)
            count = len(self.occ[len(w)])
            self.energy[offset:offset+count] = sum(values[a] for a in w)
            offset += count
        self.edges = []
        N = self.M + 1
        for mode in range(rank):
            for a, capacity in enumerate(counts):
                src, dst, coefficient = [], [], []
                for wi, w in enumerate(self.words):
                    left = capacity - w.count(a)
                    if not left:
                        continue
                    child = self.word_index[w + (a,)]
                    for j, o in enumerate(self.occ[len(w)]):
                        target = tuple(sorted(o + (mode,)))
                        src.append(self.offset[wi] + j)
                        dst.append(self.offset[child] + self.index[len(w)+1][target])
                        coefficient.append(-sqrt(left * (o.count(mode)+1)) / N)
                self.edges.append((mode, np.asarray(src, int), np.asarray(dst, int), np.asarray(coefficient)))
        self.ann = []
        self.auxoff = np.r_[0, np.cumsum([self.nwords[k] * len(self.occ[k-1]) if k else 0
                                          for k in range(self.M+1)])]
        for mode in range(rank):
            src, dst, factor = [], [], []
            number_at_k = [0] * (self.M + 1)
            for wi, w in enumerate(self.words):
                k = len(w)
                if not k:
                    continue
                base = self.auxoff[k] + number_at_k[k] * len(self.occ[k-1])
                number_at_k[k] += 1
                for j, o in enumerate(self.occ[k]):
                    count = o.count(mode)
                    if count:
                        oo = list(o)
                        oo.remove(mode)
                        src.append(self.offset[wi] + j)
                        dst.append(base + self.index[k-1][tuple(oo)])
                        factor.append(sqrt(count))
            self.ann.append((np.asarray(src, int), np.asarray(dst, int), np.asarray(factor)))
        self.rotation_groups = {}

    def block(self, c, k):
        return c[self.koff[k]:self.koff[k+1]].reshape(self.nwords[k], len(self.occ[k]))

    def hint(self, c, contact, g):
        out = np.zeros_like(c)
        for mode, src, dst, coefficient in self.edges:
            if not len(src):
                continue
            a = coefficient * (g * contact[mode])
            # Each list has distinct src/dst: last type and photon mode are fixed.
            out[dst] += a * c[src]
            out[src] += a.conj() * c[dst]
        return out

    def density(self, c):
        d = np.zeros((int(self.auxoff[-1]), self.rank), complex)
        for mode, (src, dst, factor) in enumerate(self.ann):
            d[dst, mode] = factor * c[src]
        # gamma_ij=<b_j^dagger b_i>; its eigenvectors are orbital amplitudes.
        return d.T @ d.conj()

    def field_groups(self, k, i, j):
        key = k, i, j
        if key not in self.rotation_groups:
            groups = {}
            for n, o in enumerate(self.occ[k]):
                other = tuple(x for x in o if x != i and x != j)
                pair = k - len(other)
                if pair:
                    groups.setdefault((pair, other), np.empty(pair+1, int))[o.count(i)] = n
            self.rotation_groups[key] = {n: np.stack([v for (nn, _), v in groups.items() if nn == n])
                                         for n in range(1, k+1) if any(nn == n for nn, _ in groups)}
        return self.rotation_groups[key]

    def beam(self, c, i, j, transform):
        matrices = {n: fock_two_mode(transform, n) for n in range(1, self.M+1)}
        for k in range(1, self.M+1):
            target = self.block(c, k)
            for n, indices in self.field_groups(k, i, j).items():
                # These groups partition the Fock basis; the advanced-index
                # gather is a copy, so in-place updates cannot mix groups.
                target[:, indices] = target[:, indices] @ matrices[n].T
        return c


def fock_two_mode(t, n):
    """Gamma(t) on |j,n-j>, j=0..n, including bosonic normalization."""
    out = np.zeros((n+1, n+1), complex)
    for old in range(n+1):
        b = n - old
        for new in range(n+1):
            value = 0j
            for x in range(max(0, new-b), min(old, new)+1):
                y = new - x
                value += (comb(old, x) * comb(b, y)
                          * t[0, 0]**x * t[1, 0]**(old-x)
                          * t[0, 1]**y * t[1, 1]**(b-y))
            out[new, old] = value * sqrt(factorial(new) * factorial(n-new)
                                         / (factorial(old) * factorial(b)))
    return out


def embed_coefficients(c, old, new):
    out = np.zeros(new.size, complex)
    for k in range(old.M+1):
        mapping = [new.index[k][o] for o in old.occ[k]]
        new.block(out, k)[:, mapping] = old.block(c, k)
    return out


class ResourceLimit(RuntimeError):
    pass


def commutator_bounds(g, d=.75):
    first = g * (d*d + 2*d + sqrt(2))
    second = 2*g*g*(d+1)
    return dict(AAB=first, BBA=second, strang_constant=first/24 + second/12)


class Evolution:
    def __init__(self, counts, values, modes, g=.9, rank_cap=24,
                 state_cap=2_000_000, degree=8):
        self.counts, self.values = tuple(counts), tuple(values)
        self.N = sum(counts) + 1
        self.modes, self.g, self.rank_cap, self.state_cap = modes, g, rank_cap, state_cap
        k = np.arange(1, modes+1)
        self.frequency = -2*np.cos(np.pi*k/(modes+1))
        self.contact = sqrt(2/(modes+1))*np.sin(np.pi*k/(modes+1))
        self.cache = {}
        self.layout = self.get_layout(0)
        self.c = np.array([1+0j])
        self.w = np.zeros((modes, 0), complex)
        self.degree = degree
        self.projection_error = self.polynomial_error = self.normalization_change = 0.0
        self.orthogonalization_error = 0.0
        self.maximum_rank = self.maximum_size = 0
        self.maximum_gram_error = 0.0
        self.projections = 0

    def get_layout(self, rank):
        if rank not in self.cache:
            self.cache[rank] = Layout(self.counts, self.values, rank, self.state_cap)
            while len(self.cache) > 3:
                self.cache.pop(next(iter(self.cache)))
        return self.cache[rank]

    def free(self, dt):
        self.c *= np.exp(-1j * dt * self.layout.energy)
        self.w *= np.exp(-1j * dt * self.frequency)[:, None]

    def include_contact(self):
        z = self.w.conj().T @ self.contact
        q = self.contact - self.w @ z
        # A second modified Gram-Schmidt pass pays for floating-point loss.
        correction = self.w.conj().T @ q
        z += correction
        q -= self.w @ correction
        norm = float(np.linalg.norm(q))
        if self.w.shape[1] == self.modes or norm < 1e-13:
            return z
        rank = self.w.shape[1] + 1
        if rank > self.rank_cap:
            raise ResourceLimit(f"active rank would exceed {self.rank_cap}")
        new = self.get_layout(rank)
        self.c = embed_coefficients(self.c, self.layout, new)
        self.layout = new
        self.w = np.column_stack((self.w, q/norm))
        self.maximum_rank = max(self.maximum_rank, rank)
        self.maximum_size = max(self.maximum_size, new.size)
        return np.r_[z, norm]

    def interact(self, dt, contact):
        term, out = self.c.copy(), self.c.copy()
        for j in range(1, self.degree+1):
            term = (-1j*dt/j) * self.layout.hint(term, contact, self.g)
            out += term
        norm = float(np.linalg.norm(out))
        self.normalization_change += abs(norm-1)
        self.c = out/norm
        x = self.g*abs(dt)
        # General exponential tail; normalization is accounted separately.
        self.polynomial_error += exp(x)*x**(self.degree+1)/factorial(self.degree+1)

    def orthogonalize(self):
        rank = self.w.shape[1]
        if not rank:
            return
        gram = self.w.conj().T @ self.w
        eta = float(np.linalg.norm(gram-np.eye(rank), 2))
        self.maximum_gram_error = max(self.maximum_gram_error, eta)
        if eta > 1e-14:
            values, vectors = np.linalg.eigh(gram)
            if min(values) <= 0:
                raise ResourceLimit("receiver basis lost linear independence")
            self.w = self.w @ ((vectors/np.sqrt(values)) @ vectors.conj().T)
            self.orthogonalization_error += self.layout.M*eta*exp(self.layout.M*eta)

    def compress(self, allowed):
        rank = self.w.shape[1]
        if rank <= 1 or allowed <= 0:
            return
        rho = self.layout.density(self.c)
        values, vectors = np.linalg.eigh(rho)
        if max(0.0, float(values[0])) > allowed*allowed/2:
            return
        u = vectors[:, 0].copy()
        last = rank - 1
        for i in range(last):
            a, b = u[i], u[last]
            norm = sqrt(abs(a)**2+abs(b)**2)
            if norm < 1e-15:
                continue
            t = np.array([[b, -a], [a.conjugate(), b.conjugate()]])/norm
            self.c = self.layout.beam(self.c, i, last, t)
            self.w[:, [i, last]] = self.w[:, [i, last]] @ t.conj().T
            u[i], u[last] = 0, norm
        norm = float(np.linalg.norm(self.c))
        self.normalization_change += abs(norm-1)
        self.c /= norm
        drop = 0.0
        for k in range(1, self.layout.M+1):
            indices = [j for j, o in enumerate(self.layout.occ[k]) if last in o]
            drop += float(np.sum(abs(self.layout.block(self.c, k)[:, indices])**2))
        change = sqrt(2*drop/(1+sqrt(max(0.0, 1-drop))))
        if change > allowed + 1e-9:
            raise RuntimeError("orbital projection exceeds its preassigned allowance")
        new = self.get_layout(rank-1)
        out = np.zeros(new.size, complex)
        for k in range(self.layout.M+1):
            indices = [self.layout.index[k][o] for o in new.occ[k]]
            new.block(out, k)[:] = self.layout.block(self.c, k)[:, indices]
        self.c = out/sqrt(max(1e-300, 1-drop))
        self.w = self.w[:, :-1]
        self.layout = new
        self.projection_error += change
        self.projections += 1

    def populations(self):
        return [float(np.sum(abs(self.layout.block(self.c, k))**2))
                for k in range(self.layout.M+1)]

    def run(self, duration, strang_budget=.04, compression_budget=.15,
            max_seconds=120, progress=None, steps_override=None):
        coefficient = commutator_bounds(self.g)["strang_constant"]
        steps = steps_override or max(1, ceil(sqrt(coefficient*duration**3/strang_budget)))
        dt = duration/steps
        start = time.monotonic()
        status, reason = "COMPLETE", None
        completed = 0
        self.free(dt/2)
        for step in range(steps):
            if time.monotonic()-start > max_seconds:
                status, reason = "UNRESOLVED", "per-instance wall-time budget"
                break
            try:
                contact = self.include_contact()
            except ResourceLimit as error:
                status, reason = "UNRESOLVED", str(error)
                break
            self.interact(dt, contact)
            remaining = steps-step
            unspent = max(0.0, compression_budget-self.projection_error)
            # At the hard rank ceiling spend only already-declared remaining
            # error, rather than pretending the new photon mode is absent.
            allowed = unspent if self.w.shape[1] == self.rank_cap else unspent/remaining
            self.compress(allowed)
            self.free(dt/2 if step == steps-1 else dt)
            if step % 50 == 0 or step == steps-1:
                self.orthogonalize()
            completed = step+1
            if progress and (completed % 200 == 0 or completed == steps):
                progress(dict(step=completed, steps=steps, rank=self.w.shape[1],
                              coefficients=self.layout.size, projection_error=self.projection_error,
                              elapsed_seconds=time.monotonic()-start))
        strang_error = coefficient*duration**3/steps**2
        contact_threshold_error = self.g*duration*1e-13
        error = (strang_error+self.projection_error+self.polynomial_error
                 + self.normalization_change+self.orthogonalization_error+contact_threshold_error)
        populations = self.populations()
        threshold = max(2, ceil(self.N/2))-1
        p = sum(populations[threshold:])
        lower = max(0.0, sqrt(p)-error)**2 if status == "COMPLETE" else None
        return dict(status=status, stop_reason=reason, N=self.N, receiver_modes=self.modes,
                    nonroot_type_counts=self.counts, detuning_values=self.values, g=self.g,
                    observation_time=duration, completed_physical_steps_time=completed*dt,
                    steps=steps, completed_steps=completed, dt=dt, exponential_degree=self.degree,
                    final_active_rank=self.w.shape[1], maximum_active_rank=self.maximum_rank,
                    maximum_coefficient_dimension=self.maximum_size,
                    populations=populations, macro_probability_estimate=p,
                    full_assembly_probability_estimate=populations[-1],
                    analytic_state_error_account=error,
                    analytic_probability_lower=lower,
                    analytic_full_assembly_lower=(max(0.0, sqrt(populations[-1])-error)**2
                                                  if status == "COMPLETE" else None),
                    finite_instance_gate=bool(status == "COMPLETE" and error <= .20 and lower >= .5),
                    error_components=dict(strang=strang_error, orbital_projection=self.projection_error,
                                          exponential_tail=self.polynomial_error,
                                          contact_threshold=contact_threshold_error,
                                          observed_normalization_change=self.normalization_change,
                                          polar_basis_change_upper=self.orthogonalization_error),
                    maximum_observed_gram_error=self.maximum_gram_error,
                    projections=self.projections, elapsed_seconds=time.monotonic()-start,
                    exact_history_retained=True, extra_photon_cutoff=False,
                    incomplete_population_notice=(None if status == "COMPLETE" else
                        "Populations are intermediate only, not observations at the frozen target time"),
                    roundoff_certified=False)
