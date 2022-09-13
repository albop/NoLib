include("consumption_savings_ayiagari.jl")

using ForwardDiff # autodiff
using FiniteDiff  # last resort (finitediff)


### Equilibrium function

function equilibrium(mod::typeof(model), xx, μ, y, z; diff=false, linear=false)
    
    # let's label our arguments
    s = [NoLib.LVectorLike(merge(model.m, model.s),e)  for e in model.grid[:]]
    x = NoLib.label_GArray(model.x, xx)
    y = NoLib.LVectorLike(model.y, y)
    z = NoLib.LVectorLike(model.z, z)

    p = model.p

    E = sum( 
        (s[i].y - x[i].c)*μ[i] 
        for i in 1:length(model.grid)
    )  # aggregate savings

    if !linear
        E -= p.δ*y.K   # depreciation
    end

    if !diff
        return SVector(E)
    end

    function E_x(dx)
        dx = NoLib.label_GArray(model.x, xx)
        SVector(sum( -dx[i].c*μ[i]  for i in 1:length(model.grid) ) )
    end
    E_μ = dμ->equilibrium(mod, xx, dμ, y, z; diff=false, linear=true)
    E_y = ForwardDiff.jacobian(u->equilibrium(mod, xx, μ, u, z; diff=false), y)
    E_z = ForwardDiff.jacobian(u->equilibrium(mod, xx, μ, y, u; diff=false), z)

    return (E, E_x, E_μ, E_y, E_z)
end


#### Projection function

function projection(mod::typeof(model), y, z; diff=false)
    
    # let's label our arguments
    p = model.p
    y = NoLib.LVectorLike(model.y, y)
    z = NoLib.LVectorLike(model.z, z)

    r = 1 + p.α*(1/y.K)^(1-p.α)*exp(z.z) - p.δ
    w = (1-p.α)*exp(z.z)
    
    P = SLVector((;w, r))

    if !diff
        return P
    end

    P_y = ForwardDiff.jacobian(u->projection(mod, u, z; diff=false), y)
    P_z = ForwardDiff.jacobian(u->projection(mod, y, u; diff=false), z)

    return (P, P_y, P_z)

end


#### Transition function

# Optimality conditions


y = model.y
z = model.z
p0 = projection(model, y, z)


E, E_x, E_μ, E_y, E_z = equilibrium(model, x0, μ, y, z; diff=true)
p, p_y, p_z = projection(model, y, z; diff=true)
μ1, G_μ, G_x, G_p = NoLib.G(model, μ, x0, p0; diff=true);
r, r_1, r_2, r_p0, r_p1 = NoLib.F(model, x0, x0, p0, p0; diff=true)


#####

function residual(mod::typeof(model),y)
    p = projection(model, y)
    model.grid.g1.points .= [NoLib.cover(p, e) for e in model.grid.g1.points]
    sol = NoLib.time_iteration_3(model; verbose=false, improve=false)
    x = sol.solution
    μ = NoLib.ergodic_distribution(model, x)
    res = equilibrium(model, x, μ, y)
    res
end




s0_l = [NoLib.LVectorLike(merge(model.m, model.s),e)  for e in model.grid[:]]
x0_l = NoLib.label_GArray(model.x, sol.solution)


x = sol.solution
y = LVector(K=40.0)
μ = NoLib.ergodic_distribution(model, x)



projection(model, y)
residual(model, y)

using NoLib: cover

equilibrium(model, x, μ, y; diff=true)

@time res = projection(model, y; diff=true)



resid(u::Float64) = residual(model, SLVector(K=u))

Kvec = range(15, 20;length=20)
Rvec = [resid(e) for e in Kvec]

using Plots
plot(Kvec, Rvec)

