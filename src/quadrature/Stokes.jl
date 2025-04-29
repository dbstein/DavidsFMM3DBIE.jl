
"""
Generate the combined field singular matrix: αS + βD for the Stokes problem
    pretty much a direct wrap of the fortran code
"""
function StokesCombinedFieldLayerPotential(
        source::AbstractGeometry,
        α::Float64,
        β::Float64,
        interior::Bool,
        ε::Float64=1e-8
    )::Matrix{Float64}
    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals = get_GenericContent(source)
    # wrangle things to put into Fortran code
    dpars = Float64[α, β]
    ifinout = interior ? Int32(0) : Int32(1)
    # singular quadrature matrix
    xmat = zeros(Float64, 3*npts, 3*npts)
    # call fortran
    @ccall "libfmm3dbie.so".stok_comb_vel_matgen_(
        npatches::Ref{Int32},
        norders::Ref{Int32},
        ixyzs::Ref{Int32},
        iptype::Ref{Int32},
        npts::Ref{Int32},
        srccoefs::Ref{Float64},
        srcvals::Ref{Float64},
        ε::Ref{Float64},
        dpars::Ref{Float64},
        ifinout::Ref{Int32},
        xmat::Ref{Float64}
    )::Cvoid
    return xmat
end