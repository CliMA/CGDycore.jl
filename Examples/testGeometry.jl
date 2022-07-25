using CGDycore
using LinearAlgebra

lam1 = pi 
phi1 = 0.0
lam2 = pi + 0.3 * pi
phi2 = 0.0
lam3 = pi + 0.15 * pi
phi3 = 0.3 * pi
lam = pi + 0.15 * pi
phi = 0.15 * pi

r = 1.0

P1 = CGDycore.sphere2cart(lam1,phi1,r)
P2 = CGDycore.sphere2cart(lam2,phi2,r)
P3 = CGDycore.sphere2cart(lam3,phi3,r)
P = CGDycore.sphere2cart(lam,phi,r)

@show  P1
@show  P2
@show  P3
@show  P
@show dot(cross(P1,P2),P)
@show dot(cross(P2,P3),P)
@show dot(cross(P3,P1),P)


