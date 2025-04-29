
xyz2v(x::T, y::T, z::T) where T <: Number = SVector{3, T}(x, y, z)
xyz2v(
    x::AbstractArray{T},
    y::AbstractArray{T},
    z::AbstractArray{T}
) where T <: Number = xyz2v.(x, y, z)

v2xyz(s::SVector{3, T}) where T <: Number = (s[1], s[2], s[3])
function v2xyz(
    v::AbstractArray{SVector{3, T}}
) where T <: Number
    x = T[]
    y = T[]
    z = T[]
    for _v in v
        _x, _y, _z = v2xyz(_v)
        push!(x, _x)
        push!(y, _y)
        push!(z, _z)
    end
    return x, y, z
end

Flat2SVector(
    v::AbstractVector{T}
) where T <: Number = reinterpret(SVector{3, T}, v)
SVector2Flat(
    v::AbstractVector{SVector{3, T}}
) where T <: Number = reinterpret(T, v)
