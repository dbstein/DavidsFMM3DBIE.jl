"""
Wrapper to fmm3dbie subroutine get_sphere_npat_mem_
    computes the number of patches and points
    for a sphere of radius a, centered at c0
"""
function get_sphere_npat_mem(
    a::Float64, 
    na::Int32, 
    c0::Vector{Float64}, 
    norder::Int32, 
    iptype0::Int32
)::Tuple{Int32, Int32}
    npatches = zeros(Int, 1)
    npts = zeros(Int, 1)
    @ccall "libfmm3dbie.so".get_sphere_npat_mem_(
        a::Ref{Float64},
        na::Ref{Int32},
        c0::Ref{Float64},
        norder::Ref{Int32},
        iptype0::Ref{Int32},
        npatches::Ptr{Int32},
        npts::Ptr{Int32}
    )::Cvoid
    return (npatches[1], npts[1])
end

"""
Wrapper to fmm3dbie subroutine get_sphere_npat_
"""
function get_sphere_npat(
    a::Float64,
    na::Int32,
    c0::Vector{Float64},
    norder::Int32,
    iptype0::Int32,
    npatches::Int32,
    npts::Int32
)::GenericGeometry
    gg = GenericGeometry(npatches, npts)
    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals = get_GenericContent(gg)
    # call out to fortran
    @ccall "libfmm3dbie.so".get_sphere_npat_(
        a::Ref{Float64},
        na::Ref{Int32},
        c0::Ref{Float64},
        norder::Ref{Int32},
        iptype0::Ref{Int32},
        npatches::Ref{Int32},
        npts::Ref{Int32},
        norders::Ptr{Int32},
        ixyzs::Ptr{Int32},
        iptype::Ptr{Int32},
        srccoefs::Ptr{Float64},
        srcvals::Ptr{Float64},
    )::Cvoid
    return gg
end

"""
Discretize a sphere of radius r, centered at c, with
    discretization size h and order norder.
    The discretization is done using triangular patches
    with RV nodes.
"""
function DiscretizeSphere(
    r::Float64,
    c::Tuple{Float64, Float64, Float64},
    h::Float64;
    order::Int=12,
    panel_type::Symbol=:TriangularRV,
)
    o = Int32(order)
    q = WranglePanelType(panel_type)
    c = collect(c)
    # estimate the number of patches on the rectangle
    np = ceil(Int32, ((2r / h) / order))
    # estimate the number of patches & points on the sphere
    npatches, npoints = get_sphere_npat_mem(r, np, c, o, q)
    # generate GenericGeometry object
    gg = get_sphere_npat(r, np, c, o, q, npatches, npoints)
    return PrecomputedGeometry(gg)
end

