include("../src/NoLib.jl")

cg = Main.NoLib.CartesianGrid((:a, :b), (0.0, 0.0), (1.0, 1.0), (1000, 10000))

using Main.NoLib
using LabelledArrays
using StaticArrays
using ResumableFunctions


function iter(cg::Main.NoLib.CartesianGrid{d}) where d
    # ranges = (range(e[1], e[2]; length=e[3]) for e in zip(cg.min, cg.max, cg.n))
    # ranges = tuple( (range(cg.min[i], cg.max[i];length=cg.n[i]) for i=1:d)... )
    ranges = ( range(cg.min[i], cg.max[i]; length=cg.n[i]) for i=1:length(cg.min) )
    return ( SVector(v...) for v in Iterators.product( ranges...) )  #:: Vector{Tuple{Float64, Float64}}
end

function iter3(cg::Main.NoLib.CartesianGrid{d,syms}) where d where syms
#     # ranges = (range(e[1], e[2]; length=e[3]) for e in zip(cg.min, cg.max, cg.n))
#     # ranges = tuple( (range(cg.min[i], cg.max[i];length=cg.n[i]) for i=1:d)... )
#     # ranges = ( range(cg.min[i], cg.max[i]; length=cg.n[i]) for i=1:length(cg.min) )
    ranges = ( [range(cg.min[1], cg.max[1]; length=cg.n[1])...], [range(cg.min[2], cg.max[2]; length=cg.n[2])...] )
#     tp = eltype(cg)
    tp = SLArray{Tuple{2}, Float64, 1, 2, (:a, :b)}
    return ( tp(v[1], v[2]) for v in Iterators.product( ranges...) )  #:: Vector{Tuple{Float64, Float64}}
end

@resumable function iter2(cg::Main.NoLib.CartesianGrid{d,syms}) where d where syms
    ranges = ( range(cg.min[1], cg.max[1]; length=cg.n[1]), range(cg.min[2], cg.max[2]; length=cg.n[2]) )
    tp = SLArray{Tuple{2}, Float64, 1, 2, (:a, :b)}
    # tp = eltype(cg)
    for v in Iterators.product(ranges...)
        x = tp(v...)
        @yield x
    end
end



function otest(cg::Main.NoLib.CartesianGrid{d,syms}) where d where syms
    ranges = ( range(cg.min[1], cg.max[1]; length=cg.n[1]), range(cg.min[2], cg.max[2]; length=cg.n[2]) )
    # tp = SLArray{Tuple{2}, Float64, 1, 2, (:a, :b)}
    tp = eltype(cg)
    s = 0.0
    for v in Iterators.product(ranges...)
        x = tp(v...)
        s += x.a + x.b
    end
    s
end

function mysum2(it)
    res = 0.0
    for e in it
        res += e.a+e.b
    end

    return res
end

function mysum(it)
    res = 0.0
    for e in it
        res += e[1] + e[2]
    end

    return res
end

function test(cg::Main.NoLib.CartesianGrid{d}) where d
    return mysum(iter(cg)) 
end

function test2(cg::Main.NoLib.CartesianGrid{d}) where d
    return mysum2(iter2(cg)) 
end

function test3(cg::Main.NoLib.CartesianGrid{d}) where d
    return mysum2(iter3(cg)) 
end


test2(cg)


###

using ResumableFunctions

import Base: +
import Base: length
import Base: iterate

struct GArray{T}
    grid
    data::T
end

(+)(g1::GArray, g2::GArray) = GArray(g1.grid, g1.data + g2.data)

length(g::GArray) = length(g.data)

# @resumable function iter2(cg::Main.NoLib.CartesianGrid{d,syms}) where d where syms
#     ranges = ( range(cg.min[1], cg.max[1]; length=cg.n[1]), range(cg.min[2], cg.max[2]; length=cg.n[2]) )
#     tp = SLArray{Tuple{2}, Float64, 1, 2, (:a, :b)}
#     # tp = eltype(cg)
#     for v in Iterators.product(ranges...)
#         x = tp(v...)
#         @yield x
#     end
# end

iterate(g::GArray) = iterate(g.data)
iterate(g::GArray, i) = iterate(g.data, i)


grid = nothing
data = rand(10)

x1 = GArray(grid, data)
x2 = GArray(grid, data)


f(a::Float64,b::Float64) = sin(a)+sin(b)

f.(x1, x2)



