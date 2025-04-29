"""
Single source --> Single or SIMD_WIDTH targets
"""

# computation for inverse displacement and scaled displacement vector
@inline function kernel_basics(source, target)
    @fastmath begin
        d = target .- source
        ir = 1/norm(d)
        dn = d * ir            
    end
    return ir, dn
end

# action of one stokeslet on 1/N targets
@inline function _StokesletKernel(stokeslet, ir, dn)
    scale = 0.039788735772973836 # 1/(8π)
    # scale = 1/(8π) # note this is a Float64 // required to use GPU with this kernel
    @fastmath begin
        w = dot(dn, stokeslet) * dn
        ret = (w + stokeslet) * (ir*scale)
    end
    return ret
end
@inline function StokesletKernel(source, target, stokeslet)
    ir, dn = kernel_basics(source, target)
    return _StokesletKernel(stokeslet, ir, dn)
end

# action of one stresslet on 1/N targets
@inline function _StressletKernel(stresslet, stressvec, ir, dn)
    scale = -0.238732414637843
    # scale = -3/(4π)  # note this is a Float64 // required to use GPU with this kernel
    @fastmath begin
        w1 = dot(dn, stresslet)
        w2 = dot(dn, stressvec)
        ret = dn * (w1*w2*ir*ir*scale)
    end
    return ret
end
@inline function StressletKernel(
    source, target, stresslet, stressvec)
    ir, dn = kernel_basics(source, target)
    return _StressletKernel(stresslet, stressvec, ir, dn)
end

# action of one doublet on 1/N targets
@inline function _DoubletKernel(doublet, ir, dn)
    scale = 0.039788735772973836 # 1/(8π)
    # scale = 1/(8π)  # note this is a Float64 // required to use GPU with this kernel
    @fastmath begin
        w = dot(dn, doublet) * dn
        ret = (-3w + doublet) * (ir*ir*ir*scale)
    end
    return ret
end
@inline function DoubletKernel(source, target, doublet)
    ir, dn = kernel_basics(source, target)
    return _DoubletKernel(doublet, ir, dn)
end

# Single threaded, for now
function StokesKernel!(
    velocities::AbstractArray{SVector{3, T2}},
    sources::AbstractArray{SVector{3, T1}};
    targets::Union{Nothing,AbstractArray{SVector{3, T1}}}=nothing,
    stokeslets::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
    stresslets::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
    stressvecs::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
    doublets::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
) where {T1 <: Real, T2 <: Real}
    self = targets === nothing
    targets = self ? sources : targets
    dostk = stokeslets !== nothing
    dostr = (stresslets !== nothing) && (stressvecs !== nothing)
    dodub = doublets !== nothing
    @inbounds for j in eachindex(targets, velocities)
        target = targets[j]
        velocity = zero(SVector{3, T2})
        @inbounds for i in eachindex(sources)
            if self && (i == j)
                continue
            end
            ir, dn = kernel_basics(sources[i], target)
            velocity += dostk ? _StokesletKernel(stokeslets[i], ir, dn) : zero(SVector{3, T2})
            velocity += dostr ? _StressletKernel(stresslets[i], stressvecs[i], ir, dn) : zero(SVector{3, T2})
            velocity += dodub ? _DoubletKernel(doublets[i], ir, dn) : zero(SVector{3, T2})
        end
        velocities[j] = velocity
    end
    return velocities
end
function StokesKernel(
    sources::AbstractArray{SVector{3, T1}};
    targets::Union{Nothing,AbstractArray{SVector{3, T1}}}=nothing,
    stokeslets::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
    stresslets::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
    stressvecs::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
    doublets::Union{Nothing,AbstractArray{SVector{3, T2}}}=nothing,
) where {T1 <: Real, T2 <: Real}
    self = targets === nothing
    NT = self ? length(sources) : length(targets)
    velocities = Vector{SVector{3, T2}}(undef, NT)
    return StokesKernel!(velocities, sources; targets=targets, stokeslets=stokeslets, stresslets=stresslets, stressvecs=stressvecs, doublets=doublets)
end
