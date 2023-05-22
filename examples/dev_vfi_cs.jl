using NoLib

include("models/consumption_savings.jl")


mem = NoLib.time_iteration_workspace(model)
@time sol = NoLib.time_iteration(model, mem, improve=false, interp_mode=:cubic);

nx = sol.solution


@time nx, nv = NoLib.vfi(model; improve=false, verbose=true, T=5000)

@time nx, nv = NoLib.vfi(model; improve=true, verbose=true)


s = [NoLib.enum(model.grid)...][1]
x = initial_guess(model)[1]
φ = initial_guess(model)

NoLib.F(model, s, x, φ)


using Plots


kv = [e[1] for e in model.grid[1,:]]
plot(kv,[e[1] for e in nx[1,:]]; ylims=(0, 5))
plot!(kv,[e[1] for e in nx[2,:]])
plot!(kv,kv)
# scatter!(model.calibration.s, model.calibration.x)


plot(kv,[e[1] for e in nv[1,:]]; ylims=(-2, 0))
plot!(kv,[e[1] for e in nv[2,:]])


fobj = u->-Q(model, m, s, SVector(u...), φv)

using ForwardDiff

ForwardDiff.hessian(fobj, x)