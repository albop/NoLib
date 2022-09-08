include("rbc_def.jl")


using NoLib
using NoLib: GArray
using StaticArrays
using NoLib: F, dF, F!, dF!, dF0, dF2
# using NoLib: dF2


norm(v::GArray) = maximum(u->maximum(abs, u), v.data)

function solve(model)
    T=500
     K=10
      tol_ε=1e-8
       tol_η=1e-6
        verbose=false

    N = length(model.grid)
    x0 = GArray(model.grid, [SVector(model.x) for n=1:N])
    x1 = deepcopy(x0)

    dx = deepcopy(x0)
    r0 = x0*0
    J = dF0(model, x0, x0)[2]

    local x0
    local x1

    for t=1:T
        # r0 = F(model, x0, x0)
        F!(r0, model, x0, x0)
        ε = norm(r0)
        if ε<tol_ε
            break
        end
        if verbose
            println("ϵ=$(ε)")
        end
        x1.data .= x0.data
        for k=1:K
            # r = F(model, x1, x0)
            F!(r0, model, x1, x0)
            # J = dF(model, x1, x0)
            dF!(J, model, x1, x0)
            # dx = J\r0
            for n=1:length(r0)
                dx.data[n] = J.data[n]\r0.data[n]
            end
            e = norm(dx)
            # println("e=$(e)")
            x1.data .-= dx.data
            if e<tol_η
                break
            end
        end
        x0 = x1

    end
    return x0
end

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