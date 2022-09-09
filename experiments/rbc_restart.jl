include("rbc_def.jl")

using StaticArrays

grid = model.grid

l = [SVector(model.x) for e in 1:length(grid)]

controls = GArray(
    grid,
    l
)

using NoLib: enum
using StaticArrays
#
m = SVector(grid[1][1:1]...)
s = SVector(grid[1][2:2]...)

a = controls[1]

ss  = ( (1,2), (m,s))

gen = NoLib.τ(model, ss::Tuple, a::SVector)

[gen...]


gen = NoLib.τ_fit(model, ss::Tuple, a::SVector)

[gen...]


μ = GArray(grid, ones(length(grid))/length(grid))

@time μ1 = NoLib.G(model, μ, x0)


NoLib.transition_matrix(model, x0)




sum(μ1)

[NoLib.τ_fit(model,(ss[1][1], ss[2]), a)...]