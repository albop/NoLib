module Interpolation

import LoopVectorization: VectorizationBase
import Base: getindex

getindex(A::Vector{Float64}, i::VectorizationBase.Vec{4,Int64}) = VectorizationBase.Vec{4, Float64}(A[i(1)], A[i(2)], A[i(3)], A[i(4)])


@inline function interp(ranges::Tuple{Tuple{Float64, Float64, Int64}}, values::AbstractArray{T}, x::U) where T where U 
# function interp(ranges::Tuple{Tuple{Float64, Float64, Int64}}, values, x) where T where U 

    a = (ranges[1][1])
    b = (ranges[1][2])
    n = (ranges[1][3])

    δ = (b-a)/(n-1)

    i = div( (x-a), δ)

    i = max(min(i, n-2), 0)

    λ = (x-(a + δ*i))/δ

    i_ = floor(Int, i) + 1

    v0 = values[i_]
    v1 = values[i_+1]

    (1.0 - λ)*v0 + λ*v1

end

function vecinterp_1(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    for n=1:N
        out[n] = interp(ranges, values, x[n])
    end
    return out
end


function vecinterp_2(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    out .= (u->interp(ranges, values, u)).(x)
    return out
end

using LoopVectorization

function vecinterp_3(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    @simd for n=1:N
        out[n] = interp(ranges, values, x[n])
    end
    return out
end


using LoopVectorization

function vecinterp_4(ranges, values::Vector{T}, X::Vector{<:Number}) where T
    
    K = length(X)
    out = zeros(T,K)

    a = (ranges[1][1])
    b = (ranges[1][2])
    N = (ranges[1][3])

    δ = (b-a)/(N-1)

    @turbo for n=1:K

        x = X[n]

        i = div( (x-a), δ)

        i = max(min(i, N-2), 0)

        λ = (x-(a + δ*i))/δ

        i_ = floor(Int, i) + 1

        v0 = values[i_]
        v1 = values[i_+1]
        # v0 = getindex(values, i_)
        # v1 = getindex(values, i_+1)

        out[n] = (1-λ)*v0 + λ*v1
    end
    return out
end

function vecinterp_5(ranges, values::Vector{T}, x::Vector{<:Number}) where T
    N = length(x)
    out = zeros(T,N)
    @turbo for n=1:N
        ret = interp(ranges, values, x[n])
        out[n] = ret
    end
    return out
end

end




# @code_warntype interp(ranges, values, 2.0)

# xvec = [range(-0.1, 2.1; length=1000)...]
# tvec = xvec .^ 2

# yvec = [(interp(ranges, values, x)) for x in xvec]


# using ForwardDiff

# ForwardDiff.jacobian(u->[interp(ranges, values, u[1])], [0.3])

# using Plots

# plot(xvec, tvec)
# plot!(xvec, yvec)

# using LoopVectorization: VectorizationBase


# a,b,N = ranges[1]
# x = VectorizationBase.Vec{4, Float64}(4, 3, 2, 5)
# δ = (b-a)/(N-1)

# i = div( (x-a), δ)

# i = max(min(i, N-2), 0)

# λ = (x-(a + δ*i))/δ

# i_ = floor(Int, i) + 1

# v0 = values[i_]
# v1 = values[i_+1]

# out[n] = (1-λ)*v0 + λ*v1



# import Base: getindex
# getindex(A::Vector{T}, i::VectorizationBase.Vec{4,Int64}) where T where d = VectorizationBase.Vec{4, T}(A[i(1)], A[i(2)], A[i(3)], A[i(4)])

# getindex(values, i_)


# values[i_]

# i_(1)