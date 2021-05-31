using LinearAlgebra
using HomotopyContinuation
@var x[1:3],y[1:3],z[1:3]
f0 = rand_poly([x[1],y[1],z[1]],4)
ncoef = length(support_coefficients(System([f0]))[2][1])
@var a[1:ncoef] # coeffs of quartic

f = a . monomials(4,[x[1],y[1],z[1]])

a
@var b[1:3] # coeffs of plane
surfaces = [a' * monomials([x[i],y[i],z[i]],4) for i=1:3]
planes = [b[1]*x[i]+b[2]*y[i]+b[3]*z[i]+1 for i=1:3]

D=map(i->differentiate([surfaces[i],planes[i]],[x[i],y[i],z[i]]),1:3)
detEqs = map(m->[det(m[:,1:2]),det(m[:,2:3])],D)

eqs = vcat(vcat(surfaces...),vcat(planes...),vcat(detEqs...))

S=System(eqs;parameters=a)

MR=monodromy_solve(S)
