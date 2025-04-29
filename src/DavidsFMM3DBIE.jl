module DavidsFMM3DBIE

using Bijections
using LinearAlgebra
using StaticArrays

include("geometries/GenericGeometries.jl")
export GenericGeometry, PrecomputedGeometry
export get_npts
export getX, getY, getZ, getNX, getNY, getNZ
export get_coordinates, get_normals
export getCoordinates, getNormals
export getNaiveQuad

include("geometries/Spheres.jl")
export DiscretizeSphere

include("geometries/Ellipsoids.jl")
export DiscretizeEllipsoid

include("quadrature/Stokes.jl")
export StokesCombinedFieldLayerPotential

include("kernels/Stokes.jl")
export StokesKernel!, StokesKernel

end
