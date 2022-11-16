
using NoLib

using StaticArrays
using LabelledArrays
using NoLib: SSGrid, CGrid, PGrid, GArray, DModel, LVectorLike
import NoLib: transition, arbitrage
import NoLib: ×, ⟂

# using ForwardDiff
# using FiniteDiff

using QuantEcon: rouwenhorst

model = let 

    β = 0.96
    γ = 4.0
    σ = 0.1
    ρ = 0.0
    # r = 1.025
    # y = 1.0 # available income
    
    K = 20.0
    α = 0.36
    A = 1
    δ = 0.025
    r = 1+α*(1/K)^(1-α) - δ
    w = (1-α)*K^α

    y = w

    c = 0.9*y

    λ = 0

    e = 0
    cbar = c


    m = SLVector(;w,r,e)
    s = SLVector(;y)
    x = SLVector(;c, λ)
    y = SLVector(;K)
    z = SLVector(;z=0.0)

    p = SLVector(;β, γ, σ, ρ, cbar, α, δ)

    mc = rouwenhorst(9,ρ,σ)
    
    P = mc.p
    ## decide whether this should be matrix or smatrix
    Q = [SVector(w,r,e) for e in mc.state_values] 

    N = 100

    grid = SSGrid(Q) × CGrid(((0.01,4.0,N),))
    
    name = Val(:ayiagari)

    DModel(
        (;m, s, x, y, z, p),
        grid,
        P
    )

end


function transition(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    y = exp(M.e)*M.w + (s.y-x.c)*M.r
    return SLVector( (;y) )
end


function arbitrage(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    eq = 1 - p.β*( X.c/x.c )^(-p.γ)*M.r - x.λ
    # @warn "The euler equation is satisfied only if c<w. If c=w, it can be strictly positive."
    eq2 = x.λ ⟂ s.y-x.c
    return SLVector( (;eq, eq2) )
end


sol = NoLib.time_iteration(model; verbose=false, improve=false)


# ## Plot the decision rule
# using Plots
# ga = sol.solution
# w,r,e = model.calibration.m
# cvec = [ga(2,SVector(y))[1] for y in range(0.1, 4.0; length=1000)]
# cvec = [ga(2,SVector(y))[1] for y in range(0.1, 4.0; length=1000)]
# λvec = [ga(2,SVector(y))[2] for y in range(0.1, 4.0; length=1000)]
# plot(cvec)
# plot!(λvec)




# P = NoLib.transition_matrix(model, sol.solution)



# # using Plots
# xvec = [e[1] for e in model.grid[1,:]]
# f = plot()
# for i=1:size(μ0)[1]
#     plot!(xvec,μ0[i,:])
# end
# f
# # how to plot distribution ?

# ###

# model.p.r = 1.025

# sol_2 = NoLib.time_iteration_3(model; verbose=false, improve=false)
# μ_2 = NoLib.ergodic_distribution(model, sol_2.solution)
# scatter(xvec,μ_2[1,:])
# scatter!(xvec,μ_2[2,:])
# scatter!(xvec,μ_2[3,:])
# # how to plot distribution ?

# plot([e[2] for e in model.grid[:]],μ[:])


# # sol = NoLib.time_iteration_3(model; verbose=false, improve=true)

# ###

# # Now compute the aggregate equilibrium condition




y = LVector(K=40)

function equilibrium(model, x_, μ, y; diff=false)
    p = model.calibration.p
    s = [NoLib.LVectorLike(merge(model.calibration.m, model.calibration.s),e)  for e in model.grid[:]]
    x = NoLib.label_GArray(model.calibration.x, x_)

    res = sum( μ[i]*(s[i].y-x[i].c) for i=1:length(model.grid)) - y.K*p.δ

    if diff!=true
        return [res]
    else
        res_x = dx -> sum( μ[i]*(-dx[i].c) for i=1:length(model.grid)) - y.K*p.δ
        res_μ = dμ -> sum( dμ[i]*(s[i].y-x[i].c) for i=1:length(model.grid)) - y.K*p.δ
        res_y = dy -> - dy.K*p.δ
        return res, res_x, res_μ, res_y
    end

end




## "Steady-state" values
x0 = sol.solution
μ0 = NoLib.ergodic_distribution(model, sol.solution)
p0 = SVector(model.calibration.m[1:2]...)

## Derivatives of: F

r, J_1, J_2, U, V = NoLib.F(model, x0, x0, p0, p0; diff=true);

# renormalize to get: r, I , T, U, V

r = J_1 \ r
T = J_1 \ J_2
U = J_1 \ U
V = J_1 \ V


## Derivatives of: G

ss = [NoLib.enum(model.grid)...][1]
a = x0[1]
[NoLib.τ_fit(model, ss, a, p0)...]


res = NoLib.G(model, μ0, x0, p0; diff=false)

x__ = NoLib.ravel(x0)

f = u->NoLib.ravel(NoLib.G(model, μ0, NoLib.unravel(x0, x__), p0; diff=false))

using FiniteDiff
FiniteDiff.finite_difference_jacobian(f, x__)

NoLib.G(model, μ0, x0, p0; diff=false)

using FiniteDiff


FiniteDiff.finite_difference_jacobian(
    u->NoLib.G(model, μ0, x0, p0; diff=false)
)

# renormalize to get: r, I , T, U, V

r = J_1 \ r
T = J_1 \ J_2
U = J_1 \ U
V = J_1 \ V

equilibrium(model, x0, μ0, y; diff=true)
