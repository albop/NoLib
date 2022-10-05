using NoLib: SGrid, CGrid, ×, GArray, iti
using StaticArrays

P = @SMatrix [0.9 0.1; 0.1 0.9]
Q = @SMatrix [-0.1; 0.1]

exo = SGrid( [Q[i,:] for i=1:size(Q,1)] )
endo = CGrid( ((0.1, 5.0, 100),) )
grid = exo × endo

ga = GArray(grid, [SVector(0.2, 0.1) for i=1:length(grid)])
res = [ga.data[1]]


import Base: iterate
Base.iterate(s::SGrid) = (s.points[1], 2)
Base.iterate(s::SGrid, state) = state<=length(s) ? (s.points[state], state+1) : nothing

Base.length(s::CGrid{1}) = s.ranges[1][3]
Base.iterate(s::CGrid{1}) = (s.ranges[1][1], 1)
Base.iterate(s::CGrid{1}, state) = state<length(s) ? let
    a = s.ranges[1][1]
    b = s.ranges[1][2]
    n = length(s)
    (a + state/(n-1) * (b-a), state+1)
end : nothing


function test_noalloc(exo,u=true)
    if u
        s = sum(exo)
        t = sum(s)

    else
        if typeof(exo) <: SGrid
            s = sum(exo.points)
            t = sum(s)
        elseif typeof(exo) <: CGrid{1}
            s = sum(exo.ranges[1])
            t = s
        end
    end
    if t<-100000.0
        return "HI"
    end
end

@time test_noalloc(exo, true)  # custom iterator
@time test_noalloc(exo, false) # native iterator
@time test_noalloc(endo, true)  # custom iterator
@time test_noalloc(endo, false)  # custom iterator


# @time test(exo, false) # native iterator


[endo...]
[exo...]
[grid...]


###

using NoLib: PGrid, GVector

function iterate(ga::GVector{PGrid{G1, G2, d}, T}) where d where G1<:SGrid where G2<:CGrid{1} where T
    # (SVector(ga.g1.points[1]..., ga.g2.ranges[1][1]))
    v1,s1 = iterate(ga.g1)
    v2,s2 = iterate(ga.g2)
    return (SVector(v1..., v2...), (s1, s2))
end

function iterate(ga::GVector{PGrid{G1, G2, d}, T}, state) where d where G1<:SGrid where G2<:CGrid{1} where T

    s1, s2 = state

    

    if s2 === nothing
        return nothing
    else
        if s1 === nothing
            v_1, s__1 = iterate(ga.g1)
        else
            v_1, s__1 = iterate(ga.g1, s1)
        end
        v_2, s__2 = iterate(ga.g2, s2)
    end
    # (SVector(ga.g1.points[1]..., ga.g2.ranges[1][1]))
    v1,s1 = iterate(ga.g1)
    v2,s2 = iterate(ga.g2)
    return (SVector(v1..., v2...), (s1, s2))
end



###




function test_iterate(ga, res)
    z0 = zero(eltype(ga))
    for e in ga.data
        z0 += e
    end
    res[1] = z0
    return
end

(@allocated test_iterate(ga, res)) == 0

import NoLib

function test_iterate_2(ga, res)
    z0 = zero(eltype(ga))
    for e in (ga)
        z0 += e
    end
    res[1] = z0
    return
end

(@allocated test_iterate_2(ga, res)) == 0




rr = [exo[1]]
function test_enum_exo(ga, rr)
    z0 = ga[1]
    for e in NoLib.iti(ga)
        z0 += e[2]
    end
    rr[1] = z0
    return
end
test_enum_exo(exo, rr)
@assert (@allocated test_enum_exo(exo, rr))==0


rr = [endo[1]]
function test_enum_endo(ga, rr)
    z0 = ga[1]
    for e in NoLib.iti(ga)
        z0 += e[2]
    end
    rr[1] = z0
    return
end
test_enum_endo(endo, rr)
@assert (@allocated test_enum_endo(endo, rr))==0



rr = [grid[1]]
function test_enum(ga, rr)
    z0 = ga[1]
    for e in NoLib.iti(ga)
        el = SVector(e[2][1]..., e[2][2]...)
        z0 += el
    end
    rr[1] = z0
    return
end
@time test_enum(grid, rr);