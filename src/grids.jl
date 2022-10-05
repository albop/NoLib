
abstract type AGrid{d} end
import Base: eltype, iterate

eltype(cg::AGrid{d}) where d = SVector{d, Float64}


struct CGrid{d} <: AGrid{d}
    ranges::NTuple{d, Tuple{Float64, Float64, Int64}}
end


getindex(g::CGrid{1}, i::Int) = SVector{1}(
    g.ranges[1][1] + (g.ranges[1][2]-g.ranges[1][1])*( (i-1)/(g.ranges[1][3]-1))
)

getindex(g::CGrid{1}, ::Colon) = [SVector(i) for i in range(g.ranges[1]...)]


struct SGrid{d} <: AGrid{d}
    points::Vector{SVector{d, Float64}}
end

function SGrid(Q::Matrix)
    d = size(Q,2)
    return SGrid{d}([SVector(Q[i,:]...) for i=1:size(Q,1)])
end

struct PGrid{G1, G2, d} <: AGrid{d}
    g1::G1
    g2::G2
    points::Vector{SVector{d, Float64}}
end

getindex(g::SGrid{d}, ::Colon) where d = g.points
getindex(g::SGrid{d}, i::Int) where d = g.points[i]
getindex(g::PGrid, c::CartesianIndex) = g[c[1],c[2]]

cover(m,v::SVector{d,T}) where d where T = SVector{d,T}(m...,v[length(m)+1:end]...)

PGrid(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid{SGrid{d1}, CGrid{d2}, d1+d2}(g1, g2,
    [SVector(v1...,v2...) for ((i,v1), (j,v2)) in Base.Iterators.product( iti( g1), iti(g2) ) ][:]
)

×(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid(g1,g2)

import Base: getindex
getindex(g::PGrid{G1, G2, d}, i) where G1 where G2 where d = g.points[i]
# getindex(g::PGrid{G1, G2, d}, i::Int64, j::Int64) where G1 where G2 where d = g.points[ i + length(g.g1)*(j-1)]
function getindex(g::PGrid{G1, G2, d}, i::Int64, j::Int64) where G1<:SGrid{d1} where G2<:CGrid{d2} where d where d1 where d2
    SVector{d,Float64}(g.g1[i]..., g.g2[j]...)
end

# function getit(g::PGrid{SGrid{d1}, CGrid{d2}, d}, i::Int, j::Int) where d where d1 where d2
#     SVector{d,Float64}(g.g1[i][1],g.g2[j][1])
# end


getindex(g::PGrid{G1, G2, d}, i::Int64, ::Colon) where G1 where G2 where d = g.g2[:] # TODO: should error if i out of bounds
getindex(g::PGrid{G1, G2, d}, ::Colon, i::Int64) where G1 where G2 where d = g.g1[:]


@inline to__linear_index(g::PGrid, ind::Tuple{Int64, Int64}) =  ind[1] + length(g.g1)*(ind[2]-1)


import Base: iterate
import Base: length
import Base: getindex
import Base: setindex!

iti(sg::SGrid{d}) where d = enumerate(sg.points)
# iti(cg::CGrid{d}) where d = enumerate( ( SVector(el...)    for el in Base.Iterators.product((range(r[1], r[2], r[3]) for r in cg.ranges)...) ) )

# iti2(cg::CGrid{1}) = ( (range(cg.ranges[1]...)) )
function iti(cg::CGrid{1}) 
    a = cg.ranges[1][1]
    b = cg.ranges[1][2]
    n = cg.ranges[1][3]
    ((i, SVector(a + (b-a)/(n-1)*i)) for i=0:(n-1))
end


# using ResumableFunctions

# @resumable function iti(cg::CGrid{1}) 
#     a = cg.ranges[1][1]
#     b = cg.ranges[1][2]
#     n = cg.ranges[1][3]
#     for i=0:(n-1)
#         @yield (i, SVector(a + (b-a)/(n-1)*i))
#     end
# end


function iti(pg::PGrid) 
    ( ((i,j),(v1,v2)) for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )
end


# iti(cg::CGrid{1}) = enumerate(SVector(el) for el in iti2(cg))

# function iti(pg::PGrid) 
#     ( ((i,j),(v1,v2)) for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )
# end

# enum(pg::PGrid) =  ( ( (i,j), SVector(v1...,v2...) )
#             for ((i,v1), (j,v2)) in Base.Iterators.product( enumerate(pg.g1), enumerate(pg.g2)) )


# function enum(pg::PGrid; linear_index=false) 
#     if !linear_index
#         ( ( (i,j), SVector(v1...,v2...) )
#             for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )

#     else
#         ( ( to__linear_index(pg,(i,j)), SVector(v1...,v2...) )
#             for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )        end
# end

length(pg::PGrid{G1, G2, d}) where G1 where G2 where d = length(pg.g1)*length(pg.g2)
length(sg::SGrid{d}) where d = length(sg.points)
length(cg::CGrid{d}) where d = prod(e[3] for e in cg.ranges)

### Fast iterators for 1d case

import Base: iterate

iterate(s::SGrid) = (s.points[1], 2)
iterate(s::SGrid, state) = state<=length(s) ? (s.points[state], state+1) : nothing
length(s::CGrid{1}) = s.ranges[1][3]
iterate(s::CGrid{1}) = (SVector(s.ranges[1][1]), 1)
iterate(s::CGrid{1}, state) = state<length(s) ? let
    a = s.ranges[1][1]
    b = s.ranges[1][2]
    n = length(s)
    (SVector(a + state/(n-1) * (b-a)), state+1)
end : nothing


function Base.iterate(g::PGrid)
    x,i=Base.iterate(g.g1)
    y,j=Base.iterate(g.g2)
    return (SVector(x...,y...),(i,j,y))
end
function Base.iterate(g::PGrid,state)
    i,j,y=state
    if j==length(g.g2) && i>length(g.g1)
        return nothing
    end
    if i>length(g.g1)
        i = 1
        x,i = Base.iterate(g.g1,i)
        y,j = Base.iterate(g.g2,j)
    else
        x,i = Base.iterate(g.g1,i)
    end

    return (SVector(x...,y...),(i,j,y))
end





struct EELM{I,V,O}
    value::V
    index::I
    object::O
end

# EELM(object) = EELM(Base.iterate(object)...,object)
## could be special cased for particular objects

using Base.IteratorsMD: CartesianIndices

enum(grid::PGrid) = (
    ((c[1],c[2]), grid[c]) for c in CartesianIndices((length(grid.g1), length(grid.g2)))
)

using ResumableFunctions

# enum(object::PGrid) = EELM(object)
# Base.iterate(elm::EELM{T,V,O}) where T where V where O <: PGrid = (((1,1),elm.value), elm.index)
# Base.iterate(elm::EELM{T,V,O}, state) where T where V where O <: PGrid = let 
#     it = Base.iterate(elm.object, state)
#     if it===nothing
#         return nothing
#     else
#         v,ns = it
#         return ((state[1:2],v), ns)
#     end
# end

struct C{d,d1,d2}
    grid::PGrid{SGrid{d1}, CGrid{d2},d}
    it::CartesianIndices{d}
end

eltype(cc::C) = eltype(cc.grid)
length(cc::C) = length(cc.grid)
C(grid) = C(grid, CartesianIndices((length(grid.g1), length(grid.g2))))

function iterate(cc::C)
    v, i = iterate(cc.it)
    return ( ((i[1],i[2]),cc.grid[i]) , i)
end

function iterate(cc::C, i)
    r = iterate(cc.it, i)
    if r===nothing
        return r
    end
    i = r[2]
    return ( ((i[1],i[2]),cc.grid[i]) , i)
end

enum4(grid::PGrid) = C(grid)


import Base: eachindex

eachindex(grid::PGrid{SGrid{d1}, CGrid{d2},d}) where d1 where d2 where d = ((c[1],c[2]) for c in CartesianIndices((length(grid.g1), length(grid.g2))))