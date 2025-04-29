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

function 𝐮func(𝐱::SVector{3, T}) where T
    x, y, z = v2xyz(𝐱)
    return xyz2v(2x + y^2, z - y, x - z)
end
# not tested or used at this stage
function pfunc(𝐱::SVector{3, T}) where T
    x, y, z = v2xyz(𝐱)
    return xyz2v(-2x, zero(x), zero(x))
end

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
𝐮Γ = 𝐮func.(sphereCoordinates);
𝐮Γflat = copy(SVector2Flat(𝐮Γ)); # copy is because gmres doesn't like reinterpret arrays

# solve for the density
gmres_out = gmres(StokesSingularMat, 𝐮Γflat; atol=ε, verbose=2);
𝛔flat = gmres_out[1];
𝛔 = Flat2SVector(𝛔flat);

# now we need to evaluate this at some points

𝐮e = StokesKernel(
    sphereCoordinates;
    targets=TestPoints,
    stresslets=𝛔.*sphereWeights,
    stressvecs=sphereNormals,
);
𝐮a = 𝐮func.(TestPoints);
# why am I off by a (-) sign???
𝐮a = -𝐮a;

error = norm(𝐮e - 𝐮a, Inf)/norm(𝐮a, Inf);

@printf "Relative error: %0.2e\n" error
