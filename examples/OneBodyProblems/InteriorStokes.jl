using Revise
using DavidsFMM3DBIE
using StaticArrays
using LinearAlgebra
using Printf
using Krylov

"""
Interior Dirichlet Stokes problem

    let u = (2x + y²) x̂ + (z - y) ŷ + (x - z) ẑ
        then ∇·u = 0
        and -Δu + ∇p = 0
        with p = 2x + C

This example does not use any close-evaluation, it
    is primarily a test of the singular quadratures
    and geometry objects
"""

function Ufunc(X::SVector{3, T}) where T
    x, y, z = X[1], X[2], X[3]
    return SVector{3,T}(2x + y^2, z - y, x - z)
end
# not tested or used at this stage
Pfunc(X::T) where T = -2X + one(T)

# discretize sphere
h = 0.25
panel_type = :TriangularRV
order = 12
sphere = DiscretizeSphere(2.0, (0.0, 0.0, 0.0), h; order=order, panel_type=panel_type)
Nsphere = get_npts(sphere)
# test points (near center of sphere)
TestPoints = rand(SVector{3, Float64}, 100) .- [SVector{3, Float64}(0.5, 0.5, 0.5)];

# extract some things from the sphere
sphereCoordinates = getCoordinates(sphere);
sphereNormals = getNormals(sphere);
sphereWeights = getNaiveQuad(sphere);

# generate Singular matrix
ε = 1e-8
StokesSingularMat = StokesCombinedFieldLayerPotential(sphere, 0.0, 1.0, true, ε);

# generate some boundary conditions
BoundaryU = Ufunc.(sphereCoordinates);
FlatBoundaryU = copy(reinterpret(Float64, BoundaryU));

# solve for the density
out = gmres(StokesSingularMat, FlatBoundaryU; atol=ε, verbose=2);
σ = reinterpret(SVector{3,Float64}, out[1]);

# now we need to evaluate this at some points

Uest = StokesKernel(
    sphereCoordinates;
    targets=TestPoints,
    stresslets=σ.*sphereWeights,
    stressvecs=sphereNormals,
);
Utrue = Ufunc.(TestPoints);
# why am I off by a (-) sign???
Utrue = -Utrue;

error = norm(Uest - Utrue, Inf)/norm(Utrue, Inf);

@printf "Relative error: %0.2e\n" error
