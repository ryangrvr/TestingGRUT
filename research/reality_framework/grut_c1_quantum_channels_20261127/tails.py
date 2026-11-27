"""Analytical cutoff bounds with Bernstein envelopes of radial profiles."""
import numpy as np


def cubic(profile, x):
    return np.array([profile.spline(x, k)[[0, 2]] for k in range(4)])


def bernstein(coeff, h, pad=1e-10):
    d, c, b, a = coeff
    pts = np.array([d, d+c*h/3, d+2*c*h/3+b*h*h/6, d+c*h+b*h*h/2+a*h**3/6])
    return np.array([np.min(pts, axis=0)-pad, np.max(pts, axis=0)+pad])


def multiply(a, b):
    ends = np.array([x*y for x in a for y in b])
    return np.array([np.min(ends, axis=0), np.max(ends, axis=0)])


def square(a):
    lo, hi = a
    return np.array([np.where(lo*hi <= 0, 0, np.minimum(lo**2, hi**2)), np.maximum(lo**2, hi**2)])


def potential_intervals(f2, s2m1, fs):
    return np.array([s2m1+.3*f2, s2m1+.1*f2, 3*s2m1+f2, .01*s2m1, 2*fs])


def matrix_envelope(rows):
    mids = np.mean(rows, axis=1).T
    widths = ((rows[:, 1]-rows[:, 0])/2).T
    coupled = np.zeros((len(mids), 2, 2))
    coupled[:, 0, 0], coupled[:, 1, 1] = mids[:, 0], mids[:, 2]
    coupled[:, 0, 1] = coupled[:, 1, 0] = mids[:, 4]
    radius = np.max(widths[:, :4] + widths[:, 4, None]*np.array([1, 0, 1, 0]), axis=1)
    diagonal = np.max(np.abs(mids[:, [1, 3]]), axis=1)
    op = np.maximum(np.max(np.abs(np.linalg.eigvalsh(coupled)), axis=1), diagonal) + radius
    fro_squared = np.max(np.abs(rows), axis=1)**2
    return op, np.array([1, 1, 1, 1, 2]) @ fro_squared


def edges_for(support, *profiles):
    edges = np.unique(np.concatenate([*[
        p.spline.x[(p.spline.x >= 0) & (p.spline.x < support)] for p in profiles
    ], [support]]))
    return edges[:-1], edges[1:]


def enclosure(profile, support):
    lo, hi = edges_for(support, profile)
    f, s = bernstein(cubic(profile, lo), hi-lo).transpose(1, 0, 2)
    op, fro2 = matrix_envelope(potential_intervals(square(f), square(s)-1, multiply(f, s)))
    return {
        "Vmax": np.max(op),
        "I2": 4*np.pi * np.dot(fro2, (hi**3-lo**3)/3),
        "W": np.max(hi**2*op),
        "J3": np.dot(fro2, (hi**4-lo**4)/4),
    }


def angular_sum(L, W):
    # For ell >= J+1, denominator >= 3 ell^4 and the norm ratio <= W/J^2.
    j = 100000
    ell = np.arange(L+1, j+1, dtype=float)
    r = W/(ell*(ell+1))
    if r[0] >= 1:
        raise ValueError("First omitted angular operator has no certified logarithm series")
    return np.sum(r/((1-r)*(2*ell+3)*(2*ell-1))) + W/(9*(1-W/j**2)*j**3)


def absolute(profile, support, L, Z, caps):
    inputs = enclosure(profile, support)
    m, i, w, j = (caps[k] for k in ("Vmax", "I2", "W", "J3"))
    a = 1-profile.omega**2
    r = m/(Z*Z+a)
    uv = m*i/(192*np.pi**2) / (2*Z*Z*(1-r))
    angular = Z*j/(3*np.pi)*angular_sum(L, w)
    return {"enclosure": inputs, "caps_pass": all(inputs[k] < caps[k] for k in caps),
            "frequency": uv, "angular": angular, "total": uv+angular}


def paired(first, second, support, L, Z, caps, delta):
    lo, hi = edges_for(support, first, second)
    a, b = (cubic(p, lo) for p in (first, second))
    f1, s1 = bernstein(a, hi-lo).transpose(1, 0, 2)
    f2, s2 = bernstein(b, hi-lo).transpose(1, 0, 2)
    df, ds = bernstein(a-b, hi-lo, pad=2e-10).transpose(1, 0, 2)
    changes = potential_intervals(multiply(df, f1+f2), multiply(ds, s1+s2),
                                  multiply(df, s1)+multiply(f2, ds))
    change_norm, change_fro2 = matrix_envelope(changes)
    op1, _ = matrix_envelope(potential_intervals(square(f1), square(s1)-1, multiply(f1, s1)))
    op2, _ = matrix_envelope(potential_intervals(square(f2), square(s2)-1, multiply(f2, s2)))
    env = np.maximum.reduce([op1, op2, change_norm/delta])
    w = np.max(hi**2*env)
    j = np.dot(env**2, (hi**4-lo**4)/4)
    diff_l2 = np.sqrt(4*np.pi*np.dot(change_fro2, (hi**3-lo**3)/3))
    m, b = caps["Vmax"], np.sqrt(caps["I2"])
    a = 1-max(first.omega**2, second.omega**2)
    d_omega = abs(first.omega-second.omega)
    sum_omega = abs(first.omega+second.omega)
    c = d_omega*np.sqrt(4+(sum_omega/Z)**2)
    uv = m*b/(16*np.pi**2*(1-m/(Z*Z+a))) * (diff_l2/(2*Z**2)+c*b/(3*Z**3))
    angular_prefactor = (
        delta*Z + 2/np.e*d_omega*(sum_omega*np.arctan(Z/np.sqrt(a))/np.sqrt(a) + np.log1p(Z*Z/a))
    )
    angular = 8*j/(2*np.pi)*angular_sum(L, w)*angular_prefactor
    return {"delta": delta, "W": w, "J": j, "difference_l2": diff_l2,
            "frequency": uv, "angular": angular, "total": uv+angular}
