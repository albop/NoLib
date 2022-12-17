
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

    mc = rouwenhorst(3,ρ,σ)
    
    P = mc.p
    ## decide whether this should be matrix or smatrix
    Q = [SVector(w,r,e) for e in mc.state_values] 

    N = 20

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

using NoLib: LinearOperator

function equilibrium(model, x_, μ, y_; diff=false)

    p = model.calibration.p
    s = [NoLib.LVectorLike(merge(model.calibration.m, model.calibration.s),e)  for e in model.grid[:]]
    L(u) = NoLib.label_GArray(model.calibration.x, u)
    y = NoLib.LVectorLike(model.calibration.y,y_)
    Ly(u) = NoLib.LVectorLike(model.calibration.y,y_)

    x = L(x_)


    res = sum( μ[i]*(s[i].y-x[i].c) for i=1:length(model.grid)) - y.K*p.δ

    if diff!=true
        return [res]
    else
        res_x = LinearOperator{typeof(x_), SVector{1,Float64}}(
            dx -> SVector(sum( μ[i]*(-L(dx)[i].c) for i=1:length(model.grid)))
        )
        res_μ = LinearOperator{typeof(μ), SVector{1,Float64}}(
            dμ -> SVector(sum( dμ[i]*(s[i].y-L(x)[i].c) for i=1:length(model.grid)) )
        )
        res_y = LinearOperator{SVector{1,Float64},SVector{1,Float64}}(
            dy -> SVector( -p.δ )
        )

        return (;_0=res, _x=res_x, _μ=res_μ, _y=res_y)
    end

end

using ForwardDiff

function projection(model, y_, z_; diff=false)
    p = model.calibration.p
    y = NoLib.LVectorLike(model.calibration.y, y_)
    z = NoLib.LVectorLike(model.calibration.z, z_)

    r = exp(z.z)*y.K^p.α
    w = exp(z.z)*y.K^(1-p.α)

    pp = SVector(w, r) # XXX: warning, this is order-sensitive

    if diff==false
        return pp
    end

    P_y = ForwardDiff.jacobian(u->projection(model, u, z_), y_)
    P_z = ForwardDiff.jacobian(u->projection(model, y_, u), z_)

    return (; _0=pp, _y=P_y, _z=P_z)

end


## "Steady-state" values
x0 = sol.solution   # GVector
μ0 = NoLib.ergodic_distribution(model, sol.solution)   # GDist
p0 = SVector(model.calibration.m[1:2]...)
y0 = SVector(40.0)
z0 = SVector(model.calibration.z...)

## Derivatives of: F

@time (r, J_1, J_2, U, V) = F = NoLib.F(model, x0, x0, p0, p0; diff=true);

# renormalize to get: r, I , T, U, V

r = J_1 \ r     # GVector  (GArray{G,SVector})
T = J_1 \ J_2   # Operator: GVector->GVector
U = J_1 \ U     # GMatrix (GArray{G,SMatrix})
V = J_1 \ V     # GMatrix (GArray{G,SMatrix})


## Derivatives of: G

(μ1, G_μ, G_x, G_p) = G = NoLib.G(model, μ0, x0, p0; diff=true)


# μ1: GDist ( (GArray{G,Float64}))
# G_μ: LinearOperator GDist->GDist 
# G_x: Matrix (  operates on flatten vectors  ) # todo, should be operator
# G_p: Matrix (  operates on flatten vectors  ) # todo, should probably be GMatrix

G._μ*μ0  # == μ1
G._x*x0 # same dim as - μ0
G._B*p0 # same dim as - μ0


## Derivatives of A:

(a, A_x, A_μ, A_y) = A =  equilibrium(model, x0, μ0, y0; diff=true)


# import Main.Temp: LinearOperator
# import Main.Temp: *

A._x*x0 # GVector->SVector
A._μ*μ0 # GDist->SVector
A._y*y0 # SVector->SVector

## Derivatives of P:

(p, P_y, P_z) = P = projection(model, y0, z0; diff=true)

P._y*y0
P._z*z0


# TEST
K = 100
# dy_vals = [SVector(y0*0.0001) for k=1:K]
# @time NoLib.J(F, G, A, P, dy_vals);

@time Rk, Tk, Sk = NoLib.J_forward(A, G, T, U, V, P, K);
@time Wk, Xk = NoLib.J_backward(A, G, P, Rk, K);
@time J_x, J_μ1, J_μ2 = NoLib.jacobian_matrixes(Tk, Wk, Xk, K);
