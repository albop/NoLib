
using NoLib

using StaticArrays
using LabelledArrays
using NoLib: SGrid, CGrid, PGrid, GArray
import NoLib: transition, arbitrage

model = let 

    α = 0.3
    β = 0.96
    γ = 2.0
    δ = 0.1

    k = ((1/β - (1-δ))/α)^(1/(α-1))
    i = δ*k
    z = 0.0

    p = (;α, β, γ, δ)
    
    m = SLVector( (; z))
    s = SLVector( (; k))
    x = SLVector( (; i))

    P = @SMatrix [0.9 0.1; 0.1 0.9]
    # Q = @SMatrix [-0.05; 0.05]
    Q = @SMatrix [-0.0; 0.0]

    exo = SGrid( [Q[i,:] for i=1:size(Q,1)] )
    endo = CGrid( ((0.1, 5.0, 100),) )
    grid = exo × endo


    (;m, s, x, p, P, Q, exo, endo, grid)

end

function transition(model::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    K = s.k*(1-p.δ) + x.i
    return SLVector( (;K) )
end

function arbitrage(model::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    c = exp(m.z)*s.k^p.α - x.i
    C = exp(M.z)*S.k^p.α - X.i
    r = p.β*(C/c)^(-p.γ)*(1-p.δ + p.α*exp(M.z)*S.k^(p.α-1)) - 1
    return SLVector( (;r) )
end