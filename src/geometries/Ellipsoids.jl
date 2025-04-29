"""
    get_ellipsoid_npat_mem(
        abc::Vector{Float64},
        nabc::Vector{Int32},
        c0::Vector{Float64},
        norder::Int32,
        iptype0::Int32
    )::Tuple{Int32, Int32}

Direct wrapper to FMM3DBIE subroutine `get_ellipsoid_npat_mem`
    computes the number of patches and points
    for an ellipsoid with semi-axes **abc**=[a, b, c],
    centered at **c0**, with discretization order **norder**,
    and patch type **iptype0**.

This function is not meant to be publicly called.
    Users should see DiscretizeEllipsoid instead.

Returns:\\
    - **npatches**: number of patches   \\
    - **npts**: number of points        \\

# Examples
```julia-repl
julia> abc = [0.5, 1.0, 2.0]
julia> nabc = Int32[4, 2, 8]
julia> c0 = [0.0, 0.0, 0.0]
julia> norder = Int32(8)
julia> iptype0 = Int32(1)
julia> npatches, npts = get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0)
(224, 10080)
```
"""
function get_ellipsoid_npat_mem(
    abc::Vector{Float64},
    nabc::Vector{Int32},
    c0::Vector{Float64}, 
    norder::Int32, 
    iptype0::Int32
)::Tuple{Int32, Int32}
    npatches = zeros(Int, 1)
    npts = zeros(Int, 1)
    @ccall "libfmm3dbie.so".get_ellipsoid_npat_mem_(
        abc::Ref{Float64},
        nabc::Ref{Int32},
        c0::Ref{Float64},
        norder::Ref{Int32},
        iptype0::Ref{Int32},
        npatches::Ptr{Int32},
        npts::Ptr{Int32}
    )::Cvoid
    return (npatches[1], npts[1])
end

"""
    get_ellipsoid_npat(
        abc::Vector{Float64},
        nabc::Vector{Int32},
        c0::Vector{Float64},
        norder::Int32,
        iptype0::Int32,
        npatches::Int32,
        npts::Int32
    )::GenericGeometry

Direct wrapper to FMM3DBIE subroutine `get_ellipsoid_npat`
    discretizes an ellipsoid with semi-axes **abc**=[a, b, c],
    centered at **c0**, with discretization order **norder**,
    and patch type **iptype0**.
    The number of patches and points, **npatches** and **npts**, respectively,
    is determined by a call to `get_ellipsoid_npat_mem`.

This function is not meant to be publicly called.
    Users should see DiscretizeEllipsoid instead.

# Examples
```julia-repl
julia> abc = [0.5, 1.0, 2.0]
julia> nabc = Int32[4, 2, 8]
julia> c0 = [0.0, 0.0, 0.0]
julia> norder = Int32(8)
julia> iptype0 = Int32(1)
julia> npatches, npts = get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0)
julia> gg = get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, npatches, npts)
Generic Geometry:
... # of patches:      224
... average order:     8.0
... total # of points: 10080
```
"""
function get_ellipsoid_npat(
    abc::Vector{Float64},
    nabc::Vector{Int32},
    c0::Vector{Float64},
    norder::Int32,
    iptype0::Int32,
    npatches::Int32,
    npts::Int32
)::GenericGeometry
    gg = GenericGeometry(npatches, npts)
    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals = get_GenericContent(gg)
    # call out to fortran
    @ccall "libfmm3dbie.so".get_ellipsoid_npat_(
        abc::Ref{Float64},
        nabc::Ref{Int32},
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
    DiscretizeEllipsoid(
        abc::Tuple{Float64, Float64, Float64},
        c::Tuple{Float64, Float64, Float64},
        h::Float64,
        order::Int=12,
        quadrature::Symbol=:TriangularRV,
    )::PrecomputedGeometry

Discretizes an ellipsoid with semi-axes **abc**=[a, b, c],
    centered at **c**, with approximate mesh-spacing **h**, discretization order **order**,
    and panel-type **panel\\_type**.

# Examples
```julia-repl
julia> abc = (0.5, 1.0, 2.0)
julia> c = (0.0, 0.0, 0.0)
julia> h = 0.1
julia> ell = DiscretizeEllipsoid(abc, c, h)
AbstractGeometry with type: PrecomputedGeometry
... # of patches:      56
... average order:     12.0
... total # of points: 5096
```
"""
function DiscretizeEllipsoid(
    abc::Tuple{Float64, Float64, Float64},
    c::Tuple{Float64, Float64, Float64},
    h::Float64;
    order::Int=12,
    panel_type::Symbol=:TriangularRV,
)::PrecomputedGeometry
    o = Int32(order)
    q = WranglePanelType(panel_type)
    c = collect(c)
    abc = collect(abc)
    # estimate the number of patches on the rectangle
    na = ceil(Int32, ((2abc[1] / h) / order))
    nb = ceil(Int32, ((2abc[2] / h) / order))
    nc = ceil(Int32, ((2abc[3] / h) / order))
    nabc = Int32[na, nb, nc]
    # estimate the number of patches & points on the sphere
    npatches, npoints = get_ellipsoid_npat_mem(abc, nabc, c, o, q)
    # generate GenericGeometry object
    gg = get_ellipsoid_npat(abc, nabc, c, o, q, npatches, npoints)
    return PrecomputedGeometry(gg)
end

