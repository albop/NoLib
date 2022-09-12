
using NoLib

using StaticArrays
using LabelledArrays
using NoLib: SGrid, CGrid, PGrid, GArray
import NoLib: transition, arbitrage
import NoLib: ⊗, ⟂

using QuantEcon: rouwenhorst

model = let 

    β = 0.96
    γ = 4.0
    σ = 0.1
    ρ = 0.0
    r = 1.02
    w = 1.0
    c = 0.9*w
    
    y = 0
    cbar = c


    m = SLVector(;y)
    s = SLVector(;w)
    x = S
    LVector(;c)
    p = LVector(;β, γ, σ, ρ, r, cbar)

    mc = rouwenhorst(3,ρ,σ)
    
    P = mc.p
    
    ## decide whether this should be matrix or smatrix
    Q = [SVector(e) for e in mc.state_values] 

    N = 50

    grid = SGrid(Q) ⊗ CGrid(((0.01,4.0,N),))
    
    name = Val(:consumption_savings)

    (;name, m, s, x, p, P, Q, grid)


end


function transition(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    w = exp(M.y) + (s.w-x.c)*p.r
    return SLVector( (;w) )
end

⟂(a,b) = min(a,b)
# ⟂ᶠ(a,b)

function arbitrage(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    lhs = 1 - p.β*( X.c/x.c )^(-p.γ)*p.r # >= 0
    rhs = s.w - x.c # >= 0
    res = lhs ⟂ rhs
    return SLVector( (;res) )
end



(;m,s,x,p) = model
M = m
S = s
X = x

transition(model, m, s, x, M, p)

arbitrage(model, m, s, x, M, S, X, p)




sol = NoLib.time_iteration_3(model; verbose=false, improve=false)

μ = NoLib.ergodic_distribution(model, sol.solution)

xvec = [e[1] for e in model.grid[1,:]]
plot(xvec,μ[1,:])
plot!(xvec,μ[2,:])
plot!(xvec,μ[3,:])
# how to plot distribution ?





# sol = NoLib.time_iteration_3(model; verbose=false, improve=true)




using StaticArrays

x0 = sol.solution

xval = [range(0.1, 10.0; length=100)...]

yval = [x0(2, SVector(e))[1] for e in xval]

using Plots
# plot(xval, xval)
plot(xval, xval; ylim=(0, 1.5))
plot!(xval, yval)



grid = model.grid

l = [SVector(model.x) for e in 1:length(grid)]

x0 = GArray(
    grid,
    l
)


grid = model.grid

using NoLib: enum
using StaticArrays
#
m = SVector(grid[1][1:3]...)
s = SVector(grid[1][4:4]...)

a = x0[1]

ss  = ( (1,2), (m,s))

gen = NoLib.τ(model, ss::Tuple, a::SVector)

[gen...]

NoLib.F(model, ss, x0[1], x0)