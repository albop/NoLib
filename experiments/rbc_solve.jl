include("rbc_def.jl")


using NoLib
using NoLib: GArray
using StaticArrays
using NoLib: F, dF, F!, dF!, dF0, dF2
# using NoLib: dF2


norm(v::GArray) = maximum(u->maximum(abs, u), v.data)

@code_warntype solve(model)

@time solve(model);
x = solve(model);



# @profview solve(model);

N = length(model.grid)

x0 = GArray(model.grid, [SVector(model.x) for n=1:N])
x1 = deepcopy(x0)

J = dF0(model, x0, x0)[2]

@time dF!(J, model, x1, x0)

dx = x0*0.00001

using NoLib: iti

s_ = [iti(model.grid)...][1]
x_ = x1[1]

dF2(model, s_, x_, x0, dx)

dF2(model, x1, x0, dx)


@time F!(x1, model, x0, x0)


using NoLib: ravel
using NoLib: unravel
using Plots
using LinearMaps: LinearMap
using LinearAlgebra

xx = ravel(x)

x1 = x
x0 = x
L = LinearMap(u->ravel(dF2(model, x1, x0, unravel(x0, u))), length(xx), length(xx))

dx = xx*0.0001
unravel(x0, dx)
ravel(dF2(model, x1, x0, unravel(x0, dx)))

L*(xx*0.001)

M0 = convert(Matrix, J)
M1 = convert(Matrix, L)

# plot( spy(M0.!=0), spy(M1.!=0) )

MM = M0 \ M1

maximum(abs, eigvals(MM))

################################

import NoLib: τ

s_
φ = x0

[τ(model, s_, φ)...]


S = model.grid[1]

using NoLib


g = NoLib.CGrid(((1.0, 2.0, 11),))

NoLib.trembling__hand(g, [14.000])

using NoLib: τ_fit

@time τ_fit(model, s_, x0[1])



[τ_fit(model, s_, x0[1])...]

using NoLib: GDist

using NoLib: getindex

μ = GArray(model.grid, ones(length(model.grid)))

using NoLib: G
using NoLib: setindex
G(model, μ, x0)
@code_warntype G(model, μ, x0)