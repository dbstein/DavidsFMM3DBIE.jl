"""
Solid-body rotation test problem in the unit sphere

We compute the solution to a solid-body rotation about the z-axis
    in the unit sphere.
The solution is then evaluated at the corners of a scaled unit cube,
    and compared to the exact solution
"""

using DavidsFMM3DBIE
using LinearAlgebra
using Krylov
using StaticArrays
using Printf

# how hard to make targets...
cube_const = 1/3
# tolerance for solver / evaluator
ε = 1e-8

# discretize sphere
radius = 1.0
h = 0.125
panel_type = :QuadrangularChebyshev
order = 8
sphere = DiscretizeSphere(radius, (0.0, 0.0, 0.0), h; order=order, panel_type=panel_type)
Nsphere = get_npts(sphere)

# generate the singular quad matrix
SingularMat = StokesCombinedFieldLayerPotential(sphere, 0.0, 1.0, true, ε);

# get velocity boundary condition
coordinates = getCoordinates(sphere);
X = getX(sphere);
Y = getY(sphere);
Z = getZ(sphere);
r = sqrt.(X.^2 + Y.^2);
ϕ = @. atan(Y, X); # physicists convention
uΓ = @. -sin(ϕ)*r;
vΓ = @. cos(ϕ)*r;
wΓ = zero(Z);
𝐮Γ = xyz2v(uΓ, vΓ, wΓ);
𝐮Γflat = copy(SVector2Flat(𝐮Γ)); # copy is because gmres doesn't like reinterpret arrays

# solve for density
gmres_output = gmres(SingularMat, 𝐮Γflat; atol=ε, verbose=2);
𝛔flat = gmres_output[1];
𝛔 = Flat2SVector(𝛔flat);

# construct target array
TestPoints = cube_const * [SVector{3, Float64}(i, j, k) for i in (-1.0, 1.0) for j in (-1.0, 1.0) for k in (-1.0, 1.0)]
# evaluate solution at these points
𝐮e = StokesKernel(
    getCoordinates(sphere);
    targets=TestPoints,
    stresslets=𝛔.*getNaiveQuad(sphere),
    stressvecs=getNormals(sphere),
);
# same thing, off by a (-) sign
𝐮e = -𝐮e;

# construct analytic solution
TPX, TPY, TPZ = v2xyz(TestPoints);
TPr = sqrt.(TPX.^2 + TPY.^2);
TPϕ = @. atan(TPY, TPX);
𝐮a = xyz2v(-sin.(TPϕ).*TPr, cos.(TPϕ).*TPr, zero(TPZ));

# compute error
error = norm(𝐮e - 𝐮a, Inf) / norm(𝐮a, Inf);
println("Number of points in discretization: ", Nsphere)
@printf "Error: %.2e\n" error
