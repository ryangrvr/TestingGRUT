"""Explicit review-only boundary around the sealed private kit.

No supplied source is modified. Only this facade is claimed contract-safe;
direct calls to the original module remain known unsafe and are not cleared.
The entire point sequence is materialized and validated before delegating to
appearance, nr4, certificate callbacks, or any supplied arithmetic.
"""
from fractions import Fraction


def points_snapshot(points):
    try:
        rows=tuple(dict(point) for point in points)
        weights=tuple(Fraction(row['w']) for row in rows)
    except (ValueError,TypeError,KeyError,OverflowError,ZeroDivisionError) as exc:
        raise ValueError('REFERENCE_WEIGHT_CONTRACT') from exc
    if any(w<0 for w in weights) or sum(weights)<=0:
        raise ValueError('REFERENCE_WEIGHT_CONTRACT')
    return [dict(row,w=w) for row,w in zip(rows,weights)]


class PrevalidatedKit:
    def __init__(self,supplied): self.subject=supplied

    def appearance(self,law,points,*args,**kwargs):
        return self.subject.appearance(law,points_snapshot(points),*args,**kwargs)

    def nr4(self,K,K_empty,K_S,deletions,points,*args,**kwargs):
        return self.subject.nr4(K,K_empty,K_S,deletions,points_snapshot(points),*args,**kwargs)
