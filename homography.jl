using HomotopyContinuation, LinearAlgebra, Random, GAP
# variables
@var α[1:4], β[1:4], R[1:3,1:3], t[1:2,1]
# parameters for three of the points
@var x[1:3,1:4], y[1:3,1:4]
params = vcat(vcat(x...),vcat(y...))
t = vcat(t,[1])
point_correspondences = R*x * Diagonal(α) + repeat(t,1,4)-y*Diagonal(β)
equations_SE3 = R*R'-I
planar = [det(vcat(x * Diagonal(α),ones(1,4)))]
equations = vcat(vcat(point_correspondences...),vcat([equations_SE3[i,j] for i=1:3, j=1:3 if j<=i]...),planar)
F = System(equations; parameters=params)
MR = monodromy_solve(F;permutations=true)

# use Taylor's code to investiagate the monodromy group
function ambient_Sn(MR::MonodromyResult)
    GAP.Globals.SymmetricGroup(size(permutations(MR))[1])
end
function makeGroup(MR::MonodromyResult)
    Sym = ambient_Sn(MR);
    generator_list=[];
    for i=1:size(permutations(MR))[2]
        perm = permutations(MR)[:,i];
        push!(generator_list, GAP.Globals.PermList(GAP.julia_to_gap(perm)));
    end;
    Gal = GAP.julia_to_gap(generator_list);
    GAP.Globals.Subgroup(Sym,Gal)
end
Sym = ambient_Sn(MR)
Gal = makeGroup(MR)
GAP.Globals.Order(Gal)
GAP.Globals.StructureDescription(Gal)
Aut = GAP.Globals.Centralizer(Sym,Gal);
GAP.Globals.StructureDescription(Aut)

@var H[1:3,1:3]
x_warp = H*x
correspondences = map(i->hcat(x_warp[:,i],y[:,i]),1:4)
linear_correspondences_constraints = vcat(map(m -> [det(m[1:2,:]),det(m[2:3,:])],correspondences)...)
S = System(vcat(linear_correspondences_constraints, [det(H'*H-I)]);parameters=params)

# for this formulation, it helps to fabricate a start pair
# generate calibrated homography matrix w/ normal N
function randomHomographyMatrix(N::Array{Float64})
    A = rand(3,3);
    R = (I+A-A')*inv(I-A+A');
    R + rand(3,1) * N
end
x_init = rand(3,4)
N_init = rand(1,3)
H_init = randomHomographyMatrix(N_init)
y_init = H_init * x_init
solution_init = reshape(H_init,9)
parameter_init = vcat(reshape(x_init,12),reshape(y_init,12))
sign_symmetry(h) = [h,-h]
M=monodromy_solve(S,solution_init,parameter_init;group_action=sign_symmetry)


## two homographies
@var H1[1:3,1:3]
@var H2[1:3,1:3]

# "obvious" constraints on homographies
W1 = H1'*H1-I
W2 = H2'*H2-I
F=[det(W1),det(W2),det(W1+W2),det(W1-W2)]

# parametric constraints
@var x[1:3,1:3], cx[1:2,1]
@var y[1:3,1:3], cy[1:2,1]
@var z[1:3,1:3], cz[1:2,1]
params = vcat(reshape(x,9),reshape(y,9),reshape(z,9),cx,cy,cz...)
x = hcat(x, x[:,1:2]*cx)
y = hcat(y, y[:,1:2]*cy)
z = hcat(z, z[:,1:2]*cz)
correspondence12 = vcat(map(i->cross(H1*x[:,i],y[:,i])[1:2], 1:4)...)[1:7]
correspondence13 = vcat(map(i->cross(H2*x[:,i],z[:,i])[1:2], 1:4)...)[1:7]
S=System(vcat(F,correspondence12,correspondence13);parameters=params)

# fabricate data
x_init = rand(3,3)
cx_init = rand(2,1)
x_init = hcat(x_init, x_init[:,1:2] * cx_init)
N_init = rand(1,3)
H1_init = randomHomographyMatrix(N_init)
y_init = H1_init * x_init
ny_init = nullspace(hcat(y_init[:,1:2],y_init[:,4]))
cy_init = (-1/ny_init[3]) * ny_init[1:2]
H2_init = randomHomographyMatrix(N_init)
z_init = H2_init * x_init
nz_init = nullspace(hcat(z_init[:,1:2],z_init[:,4]))
cz_init = (-1/nz_init[3]) * nz_init[1:2]
solution_init = vcat(reshape(H1_init,9),reshape(H2_init,9))
parameter_init = vcat(reshape(x_init[:,1:3],9),reshape(y_init[:,1:3],9),reshape(z_init[:,1:3],9),reshape(cx_init,2),cy_init,cz_init...)

# no more twisted pair, but there is a reflective symmetry between each pair of views
sign_symmetries(h) = [h,-h,vcat(h[1:9],-h[10:18]),vcat(-h[1:9],h[10:18])]
M=monodromy_solve(S,solution_init,parameter_init;group_action=sign_symmetries)

# blackbox solve?
# in view 1j (j=2,3):
# point-correspondences 1, 2, & 3 give two independent linear equations
# point 4 gives an additional linear equation -- 7 total
correspondence12_init = vcat(map(i->cross(H1*x_init[:,i],y_init[:,i])[1:2], 1:4)...)[1:7]
correspondence13_init = vcat(map(i->cross(H2*x_init[:,i],z_init[:,i])[1:2], 1:4)...)[1:7]
S_specialized=System(vcat(F,correspondence12_init,correspondence13_init))
result_blackbox = solve(S_specialized;start_system = :total_degree)

# separate into non-singular and singular
good_endpoints = solutions(result_blackbox)
finite_endpoints = solution.(filter(is_success,path_results(result_blackbox)))
bad_endpoints = filter(x->~(x in good_endpoints),finite_endpoints)
maximum(norm.(S_specialized.(bad_endpoints)))

# two octic constraints filter out some of the bad_endpoints
R1 = det([
    [W1[1,1] -2*W1[1,3] W1[3,3] 0]
    [0 W1[1,1] -2*W1[1,3] W1[3,3]]
    [W2[1,1] -2*W2[1,3] W2[3,3] 0]
    [0 W2[1,1] -2*W2[1,3] W2[3,3]]
])
R2 = det([
    [W1[2,2] -2*W1[2,3] W1[3,3] 0]
    [0 W1[2,2] -2*W1[2,3] W1[3,3]]
    [W2[2,2] -2*W2[2,3] W2[3,3] 0]
    [0 W2[2,2] -2*W2[2,3] W2[3,3]]
])
S_octics = System([R1,R2])
maximum(norm.(S_octics.(good_endpoints)))

# the rest of bad_endpoints have singular homography matrices
D = System([det(H1),det(H2)])
minimum(norm.(D.(good_endpoints)))
filter(x->(norm(S_octics(x))<1e-5) & (norm(D(x))>1e-5),bad_endpoints)
filter(x->(norm(S_octics(x))<1e-5) & (norm(D(x))>1e-5),good_endpoints)
