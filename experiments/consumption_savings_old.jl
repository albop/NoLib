
using NoLib

using StaticArrays
using LabelledArrays
using NoLib: SGrid, CGrid, PGrid, GArray
import NoLib: transition, arbitrage
import NoLib: ⊗, ⟂

using QuantEcon: rouwenhorst

model = let 

    β = 0.99
    B = 1e-10
    a_max = 300.
    a = 1.0
    i = a

    ρ = 0.95
    σ = 0.06^2

    # definition or r,w
    α = 0.36
    K = 41.
    A = 1
    δ = 0.025
    r = α*(1/K)^(1-α) - δ
    w = (1-α)

    e = 1

    m = SLVector(;r, w, e)
    s = SLVector(;a)
    x = SLVector(;i)
    p = SLVector(;β, B, a_max)

    mc = rouwenhorst(3,ρ,σ)
    
    P = mc.p
    ## decide whether this should be matrix or smatrix
    Q = [SVector(r,w,e) for e in mc.state_values] 

    N = 30

    grid = SGrid(Q) ⊗ CGrid(((a,B,N),))
    
    name = Val(:consumption_savings)

    (;name, m, s, x, p, P, Q, grid)


end


function transition(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    a = x.i
    return SLVector( (;a) )
end

⟂(a,b) = min(a,b)
# ⟂ᶠ(a,b)

function arbitrage(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    c = (1+m.r)*s.a +  m.w*exp(m.e) - x.i
    C = (1+M.r)*S.a +  M.w*exp(M.e) - X.i
    lhs = 1-p.β*(1+M.r)*c/C # >=0
    rhs = x.i - (-p.B)                   # >= 0
    res = lhs ⟂ rhs
    return SLVector( (;res) )
end



(;m,s,x,p) = model
M = m
S = s
X = x

transition(model, m, s, x, M, p)

arbitrage(model, m, s, x, M, S, X, p)




sol = NoLib.time_iteration_3(model; verbose=true, improve=false)

using StaticArrays

x0 = sol.solution

xval = [range(1.0, 10.0; length=10)...]

yval = [x0(2, SVector(e))[1] for e in xval]

using Plots
# plot(xval, xval)
plot(xval, yval)



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