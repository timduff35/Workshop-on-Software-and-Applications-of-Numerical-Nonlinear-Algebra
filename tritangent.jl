using LinearAlgebra
using HomotopyContinuation
ncoef =binomial(3+4,3)
@var x[1:3],y[1:3],z[1:3] # coordinates of points
@var a[1:ncoef] # coefficients of quartic
@var b[1:3] # coefficients of equations of plane
surfaces = [a' * monomials([x[i],y[i],z[i]],4) for i=1:3]
planes = [b[1]*x[i]+b[2]*y[i]+b[3]*z[i]+1 for i=1:3]
D=map(i->differentiate([surfaces[i],planes[i]],[x[i],y[i],z[i]]),1:3)
detEqs = map(m->[det(m[:,1:2]),det(m[:,2:3])],D)
eqs = vcat(vcat(surfaces...),vcat(planes...),vcat(detEqs...))
S=System(eqs;parameters=a)
MR=monodromy_solve(S)
