
abstract type AbstractGeometry end
get_npatches(ag::AbstractGeometry) = get_npatches(get_GenericGeometry(ag))
get_npts(ag::AbstractGeometry) = get_npts(get_GenericGeometry(ag))
get_norders(ag::AbstractGeometry) = get_norders(get_GenericGeometry(ag))
get_ixyzs(ag::AbstractGeometry) = get_ixyzs(get_GenericGeometry(ag))
get_iptype(ag::AbstractGeometry) = get_iptype(get_GenericGeometry(ag))
get_srccoefs(ag::AbstractGeometry) = get_srccoefs(get_GenericGeometry(ag))
get_srcvals(ag::AbstractGeometry) = get_srcvals(get_GenericGeometry(ag))
get_GenericContent(ag::AbstractGeometry) = get_GenericContent(get_GenericGeometry(ag))
Base.print(ag::AbstractGeometry) = show(ag)
Base.display(ag::AbstractGeometry) = show(ag)
function Base.show(ag::AbstractGeometry)
    println("AbstractGeometry with type: ", typeof(ag))
    _show(get_GenericGeometry(ag))
end
function getNaiveQuad(
    ag::AbstractGeometry
)::Vector{Float64}
    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals = get_GenericContent(ag)
    qwts = zeros(Float64, npts)
    @ccall "libfmm3dbie.so".get_qwts_(
        npatches::Ref{Int32},
        norders::Ptr{Int32},
        ixyzs::Ptr{Int32},
        iptype::Ptr{Int32},
        npts::Ref{Int32},
        srcvals::Ptr{Float64},
        qwts::Ptr{Float64}
    )::Cvoid
    return qwts
end
function GetPatchId(
    ag::AbstractGeometry
)::Tuple{Vector{Int32}, Matrix{Float64}}
    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals = get_GenericContent(ag)
    ipatch_id = zeros(Int32, npts)
    uvs_pts = zeros(Float64, 2, npts)
    @ccall "libfmm3dbie.so".get_patch_id_uvs_(
        npatches::Ref{Int32},
        norders::Ref{Int32},
        ixyzs::Ref{Int32},
        iptype::Ref{Int32},
        npts::Ref{Int32},
        ipatch_id::Ref{Int32},
        uvs_pts::Ref{Float64}
    )::Cvoid
    return ipatch_id, uvs_pts
end

"""
Generic Geometry

This is really a pass-through type, meant to make working with
FMM3DBIE easier.  Users shoul use PrecomputedGeometry, instead
"""
struct GenericGeometry <: AbstractGeometry
    npatches::Int32
    npts::Int32
    norders::Vector{Int32}
    ixyzs::Vector{Int32}
    iptype::Vector{Int32}
    srccoefs::Matrix{Float64}
    srcvals::Matrix{Float64}
end
function GenericGeometry(npatches, npts)
    norders = zeros(Int32, npatches)
    ixyzs = zeros(Int32, npatches+1)
    iptype = zeros(Int32, npatches)
    srccoefs = zeros(Float64, 9, npts)
    srcvals = zeros(Float64, 12, npts)
    return GenericGeometry(npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)
end
get_GenericGeometry(gg::GenericGeometry) = gg
get_npatches(gg::GenericGeometry) = gg.npatches
get_npts(gg::GenericGeometry) = gg.npts
get_norders(gg::GenericGeometry) = gg.norders
get_ixyzs(gg::GenericGeometry) = gg.ixyzs
get_iptype(gg::GenericGeometry) = gg.iptype
get_srccoefs(gg::GenericGeometry) = gg.srccoefs
get_srcvals(gg::GenericGeometry) = gg.srcvals
get_GenericContent(gg::GenericGeometry) = (get_npatches(gg), get_npts(gg), get_norders(gg), get_ixyzs(gg), get_iptype(gg), get_srccoefs(gg), get_srcvals(gg))
function Base.show(gg::GenericGeometry)
    println("Generic Geometry:")
    _show(gg)
end
function _show(gg::GenericGeometry)
    qtype = WranglePanelType(get_iptype(gg)[1])
    println("... # of patches:      ", get_npatches(gg))
    println("... patch 1 type:      ", string(qtype))
    println("... average order:     ", sum(get_norders(gg))/get_npatches(gg))
    println("... total # of points: ", get_npts(gg))
end

"""
Precomputed Geometries for FMM3DBIE
The fields:
    npatches
    norders
    ixyzs
    iptype
    npts
    srccoefs
    srcvals
Are direct copy-overs from how Fortran handles things.
I've also precomputed naive quads, ipatch_ids, and uvs_pts
Perhaps more duplicatively, I've made a bunch of copies
    for the points, normals, coords. This makes life easier
    at this development stage but these should perhaps be removed
    for brevity. As a note though, even for a million point geometry,
    these only add 96mb of memory, so...
"""
struct PrecomputedGeometry <: AbstractGeometry
    # generic geometry
    gg::GenericGeometry
    # precomputes
    naive_quad::Vector{Float64}
    ipatch_id::Vector{Int32}
    uvs_pts::Matrix{Float64}
    # these are duplicative in memory, should decide what rep to use
    # but for now will keep a few around to see what works best
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    coordinates::Matrix{Float64}
    Coordinates::Vector{SVector{3, Float64}}
    NX::Vector{Float64}
    NY::Vector{Float64}
    NZ::Vector{Float64}
    normals::Matrix{Float64}
    Normals::Vector{SVector{3, Float64}}
end
function PrecomputedGeometry(
        gg::GenericGeometry,
    )
    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals = get_GenericContent(gg)
    X = copy(srcvals[1,:])
    Y = copy(srcvals[2,:])
    Z = copy(srcvals[3,:])
    NX = copy(srcvals[10,:])
    NY = copy(srcvals[11,:])
    NZ = copy(srcvals[12,:])
    coordinates = [X'; Y'; Z']
    Coordinates = [SVector{3, Float64}(X[i], Y[i], Z[i]) for i in 1:npts]
    normals = [NX'; NY'; NZ']
    Normals = [SVector{3, Float64}(NX[i], NY[i], NZ[i]) for i in 1:npts]
    naive_quad = getNaiveQuad(gg)
    ipatch_id, uvs_pts = GetPatchId(gg)
    return PrecomputedGeometry(
                gg,
                naive_quad, ipatch_id, uvs_pts,
                X, Y, Z, coordinates, Coordinates,
                NX, NY, NZ, normals, Normals
            )
end
get_GenericGeometry(pg::PrecomputedGeometry) = pg.gg
getX(pg::PrecomputedGeometry) = pg.X
getY(pg::PrecomputedGeometry) = pg.Y
getZ(pg::PrecomputedGeometry) = pg.Z
getNX(pg::PrecomputedGeometry) = pg.NX
getNY(pg::PrecomputedGeometry) = pg.NY
getNZ(pg::PrecomputedGeometry) = pg.NZ
get_coordinates(pg::PrecomputedGeometry) = pg.coordinates
getCoordinates(pg::PrecomputedGeometry) = pg.Coordinates
get_normals(pg::PrecomputedGeometry) = pg.normals
getNormals(pg::PrecomputedGeometry) = pg.Normals
getNaiveQuad(pg::PrecomputedGeometry) = pg.naive_quad
getIPatchID(pg::PrecomputedGeometry) = pg.ipatch_id
getUVSpts(pg::PrecomputedGeometry) = pg.uvs_pts

# Used for transforming between Symbols and Int32s for Julia vs Fortran Code
const PanelTypeWrangler = Bijection{Symbol, Int32}()
PanelTypeWrangler[:TriangularRV] = 1
PanelTypeWrangler[:QuadrangularGL] = 11
PanelTypeWrangler[:QuadrangularChebyshev] = 12
WranglePanelType(pSymbol::Symbol) = PanelTypeWrangler[pSymbol]
WranglePanelType(pInt32::Int32) = PanelTypeWrangler(pInt32)
