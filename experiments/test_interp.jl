
using NoLib.Interpolation: vecinterp_1, vecinterp_2, vecinterp_3, vecinterp_4, vecinterp_5

ranges = ( (0.0,2.0,10), )

N = 10000000
xvec = [range(0,1, N)...]
values = [s^2 for s in range(ranges[1]...)]



using BenchmarkTools

out_1 = vecinterp_1(ranges, values, xvec)
out_2 = vecinterp_2(ranges, values, xvec)
out_3 = vecinterp_3(ranges, values, xvec)
out_4 = vecinterp_4(ranges, values, xvec)
out_5 = vecinterp_5(ranges, values, xvec)


@benchmark out_1 = vecinterp_1(ranges, values, xvec)
@benchmark out_2 = vecinterp_2(ranges, values, xvec)
@benchmark out_3 = vecinterp_3(ranges, values, xvec)
@benchmark out_4 = vecinterp_4(ranges, values, xvec)
@benchmark out_5 = vecinterp_5(ranges, values, xvec)


using LoopVectorization

function interp(N,v,x)
    i = floor(Int64,x*N)
    i = max(min(i, N),1)
    return exp(getindex(v,i))
end

function example(N=10, K=100)
    v = rand(N)
    x = rand(K)
    out = zeros(K)
    @tturbo for i=1:K
        ret = interp(N, v, getindex(x,i))
        out[i] = ret
    end
    return out
end

example()

import LoopVectorization: VectorizationBase

b= VectorizationBase.Vec{4,Int64}(6,7,3,6)

import Base: getindex

getindex(A::Vector{Float64}, i::VectorizationBase.Vec{4,Int64}) = VectorizationBase.Vec{4, Float64}(A[i(1)], A[i(2)], A[i(3)], A[i(4)])



vv = rand(100)

getindex(vv, b)
vv[b]

vv[b]
example()

function example_inlined(N=10, K=100)
    vvec = rand(N)
    xvec = rand(K)
    out = zeros(K)
    @turbo for n=1:K
        x = xvec[n]
        i = floor(Int64,x*N)
        i = max(min(i, N),1)
        ret = exp(vvec[i])
        out[n] = ret
    end
    return out
end


example_inlined()