
abstract type AGrid{d} end

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

cover(m,v::SVector{d,T}) where d where T = SVector{d,T}(m...,v[length(m)+1:end]...)

PGrid(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid{SGrid{d1}, CGrid{d2}, d1+d2}(g1, g2,
    [SVector(v1...,v2...) for ((i,v1), (j,v2)) in Base.Iterators.product( iti( g1), iti(g2) ) ][:]
)

×(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid(g1,g2)

import Base: getindex
getindex(g::PGrid{G1, G2, d}, i) where G1 where G2 where d = g.points[i]
getindex(g::PGrid{G1, G2, d}, i::Int64, j::Int64) where G1 where G2 where d = g.points[ i + length(g.g1)*(j-1)]
getindex(g::PGrid{G1, G2, d}, i::Int64, ::Colon) where G1 where G2 where d = g.g2[:] # TODO: should error if i out of bounds
getindex(g::PGrid{G1, G2, d}, ::Colon, i::Int64) where G1 where G2 where d = g.g1[:]


@inline to__linear_index(g::PGrid, ind::Tuple{Int64, Int64}) =  ind[1] + length(g.g1)*(ind[2]-1)


import Base: iterate
import Base: length
import Base: getindex
import Base: setindex!

iti(sg::SGrid{d}) where d = enumerate(p for p in sg.points)
iti(cg::CGrid{d}) where d = enumerate(  SVector(el...)    for el in Base.Iterators.product((range(r[1], r[2], r[3]) for r in cg.ranges)...) ) 
iti(pg::PGrid) = ( ((i,j),(v1,v2)) for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )

function enum(pg::PGrid; linear_index=false) 
    if !linear_index
        ( ( (i,j), SVector(v1...,v2...) )
            for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )

    else
        ( ( to__linear_index(pg,(i,j)), SVector(v1...,v2...) )
            for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )        end
end

length(pg::PGrid{G1, G2, d}) where G1 where G2 where d = length(pg.g1)*length(pg.g2)
length(sg::SGrid{d}) where d = length(sg.points)
length(cg::CGrid{d}) where d = prod(e[3] for e in cg.ranges)

