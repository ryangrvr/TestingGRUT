from fractions import Fraction as F

jp,jm=F(1,64),F(1,96)
b=2*(jp+jm)
a=F(13,4)-F(5,2)
q=b/a
c=b*b/a/(1-q*q)
assert b==F(5,96)
assert a==F(3,4)
assert q==F(5,72)
assert c==F(75,20636)
assert (jp*jp+jm*jm)/F(5,4)==F(13,46080)
assert 2*jp*jm/F(5,4)==F(1,3840)
print('six exact rational identities pass; no source dynamics or operating probability evaluated')
