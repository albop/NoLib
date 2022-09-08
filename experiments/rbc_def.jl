
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
    Q = @SMatrix [-0.05; 0.05]

    exo = SGrid( [Q[i,:] for i=1:size(Q,1)] )
    endo = CGrid( ((2.0, 5.0, 100),) )
    grid = PGrid(exo, endo)


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



# s0 = [NoLib.iti(model.grid)...][1]
# x0 = SVector(model.x)
# xvec = [x0 for i=1:length(model.grid)]

# xv = GArray(model.grid, xvec)




# using BenchmarkTools



# @time NoLib.F0(model, xv, xv)
# @time NoLib.F(model, xv, xv)


# @time NoLib.F0(model, s0, x0, xv)
# @time NoLib.F(model, s0, x0, xv)


# @code_warntype NoLib.F(model, s0, x0, xv);

# @time NoLib.F0(model, xv, xv)
# @time NoLib.F(model, xv, xv)

# @benchmark NoLib.dF(model, xv, xv)
# @benchmark NoLib.dF0(model, xv, xv)

# # @time NoLib.F(model, s0, xv)

# g =  NoLib.F(model, xv, xv)
# dg =  NoLib.dF(model, xv, xv)
# gdg =  NoLib.dF(model, xv, xv)


# norm(v) = maximum(u->maximum(abs, u), v.data)
# function time_iteration(model, x0; maxit=10)
#     x1 = x0
#     for it=1:maxit
#         for k=1:5
#             dr = NoLib.F(model, x0, x1)
#             A = NoLib.dF(model, x0, x1)
#             ε = norm(dr)
#             δ = A\dr
#             η = norm(δ)
#             println("ε $ε | η $η")

#             x0 = x0 - δ
#         end
#         x1 = x0
#     end
# end