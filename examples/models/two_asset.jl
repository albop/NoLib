using NoLib
const NL=NoLib

using StaticArrays
using LabelledArrays
using NoLib: SSGrid, CGrid, PGrid, GArray, DModel
import NoLib: ×, ⟂
import NoLib: transition, arbitrage, recalibrate, initial_guess, projection, equilibrium


using QuantEcon: rouwenhorst

model = let 

    ### Parameters
    β  = 0.962  # Discount factor
    θ  = 0.749  # Disutility of labour
    σ  = 2.0    # Inverse of eis
    ν  = 1.0    # Inverse of Frisch elasticity
    α  = 0.36
    χ1 = 2.589
    χ2 = 2.0
    χ0 = 0.25
    ω  = 0.005
    δ  = 0.02
    κ  = 0.1
    ϕ  = 1.5
    ϕy = 0.0
    τ  = 0.1
    r  = 0.05

    Bh = 1.04 # Exogenous liquid assets
    Bg = 2.8  # Government debt 
    G  = 0.2  # Government spending 
    K  = 10   # Capital
    total_wealth = 14

    λ  = 0
    μa = 0
    μb = 0
    I  = δ * K
    ra = r
    rb = r - ω 

    ### Variables
    m = SLVector(;w,rb,ra,e) # Prices
    s = SLVector(;y)  # States
    x = SLVector(;c,n,a,b,λ,μa,μb,port_costs)  # Controls
    y = SLVector(;Bh,Bg,G,K,total_wealth)  # Exogenous
    z = SLVector(;z=0)

    p = SLVector(;β,γ,ν,θ,ρE,σE,χ0,χ1,χ2,a_min,b_min)  #### TO DO!!

    ### Grids
    nE = 3
    ρE = 0.966
    σE = 0.92 
    mc = rouwenhorst(nE,ρE,σE)
    
    P = convert(SMatrix{nE, nE}, mc.p)
    Q = SVector( (SVector(w,r,e) for e in mc.state_values)... )

    a_min  = 0.0
    a_max  = 400.0
    nA     = 70
    b_min  = 0.0
    b_max  = 50.0
    nB     = 50
    a_grid = SSGrid(Q) × CGrid(((a_min, a_max, nA),))
    b_grid = SSGrid(Q) × CGrid(((b_min, b_max, nB),))

    name = Val(:two_asset)

    DModel(
        (;m, s, x, y, z, p),
        (;a_grid, b_grid),
        P
    )

end


function transition(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    y = exp(M.e)*M.w + (s.y-x.c)*(1+M.r)
    return SLVector( (;y) )
end

function Psi1(m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    return sign(X.A - (1.0+m.ra)*x.a) * p.χ1 * abs( (X.A - (1.0+m.ra)*x.a)/((1.0+m.ra)*x.a - p.χ0) )^(p.χ2-1.0)
end

function Psi2(m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    return -(1.0+m.ra)*(Psi1(m,s,x,M,S,X,p) + p.χ1*(p.χ2-1.0) / p.χ2*(abs(X.A-(1.0+m.ra)*x.a)/((1.0+m.ra)*x.a+p.χ0))^p.χ2)
end

function arbitrage(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    eq  = 1 - p.β*( ((X.c^(-p.γ))*(1+M.rb)) / ((x.c^(-p.γ)) - x.μb) )
    eq2 = - p.θ * (x.n^p.ν) + (x.c^(-p.γ)) * (1.0-p.τ) * m.w * m.e
    eq3 = 1 - p.β*( ((X.c^(-p.γ))*(1+M.ra-Psi2(m,s,x,M,S,X,p)))/((x.c^(-p.γ))*(1.0+Psi1(m,s,x,M,S,X,p))+x.μa) )
    eq4 = x.μb ⟂ x.b-p.b_min
    eq5 = x.μa ⟂ x.a-p.a_min
    return SLVector( (;eq, eq2, eq3, eq4) )
end

# function equilibrium(model, s::SLArray, x::SLArray, y::SLArray, p)
#     res1 = s.y - x.c - g - i - x.port_cost - psip
#     res2 = x.a + x.b - p - y.Bg
#     res3 = x.n - y.N
#     SLVector((;res1, res2, res3))
# end

# function projection(model, y::SLArray, z::SLArray, p)
#     r = p.α*exp(z.z)*(1/y.K)^(1-p.α) - p.δ
#     w = (1-p.α)*exp(z.z)*y.K^(p.α)
#     return SLVector((;w, r)) # XXX: warning, this is order-sensitive
    
# end

# function initial_guess(model, m::SLArray, s::SLArray, p)
#     c = s.y*0.9
#     a = 
#     b = 
#     return SLVector(;c)
# end

# function reward(model::typeof(model), s::SLArray, x::SLArray, p)
#     c = x.c
#     n = x.n
#     return (c^(1.0-p.γ) / (1.0-p.γ)) - p.θ*(n^(1.0+p.ν) / (1.0+p.ν))
# end


