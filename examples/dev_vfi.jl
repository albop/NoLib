using NoLib

include("models/neoclassical.jl")



φ = NoLib.initial_guess(model)
x0 = NoLib.initial_guess(model).data


φv = NoLib.GArray(model.grid, [1.0 for i=1:length(model.grid)])

# φv = NoLib.DFun(model, v0)

s0 = [NoLib.enum(model.grid)...][5]

φv(s0)

function Q(model, s, x, φv)

    p = model.calibration.p
    r = reward(model, s, x)
    contv = sum( w*φv(S)  for (w,S) in NoLib.τ(model, s, x))
    return r + p.β*contv

end


function X(model, s)

    
    p = model.calibration.p

    y = exp(s[2][1])*s[2][2]^p.α 
    ([0], [y*0.95])

end

x = x0[5]

m = SVector(model.calibration.m...)
s = SVector(model.calibration.s...)
s = [NoLib.enum(model.grid)...][5]
x = SVector(model.calibration.x...)

Q(model, m, s, x, φv)



using Optim

xx0 = SVector(x...)

NoLib.LVectorLike(x,xx0)

fobj = u->-Q(model, m, s, SVector(u...), φv)

xv0 = Vector(xx0)
@time res = Optim.optimize(fobj, xv0)

function step(model, x0, φv)

    nx = deepcopy(x0)
    nv = deepcopy(φv)

    for (n,(s,x)) in enumerate(zip(enum( model.grid ),x0))

        lb, ub = X(model, s)

        res = Optim.optimize(
            u->-Q(model, s, SVector(u...), φv),
            lb,
            ub,
            Vector(x);
            autodiff=:forward
        )

        nx[n] = res.minimizer
        nv[n] = -res.minimum
    end

    return GVector(model.grid, nx), nv

end


@time nx,nv = step(model, x0, φv)

@time nx,nv = step(model, nx, nv)

function solve(model, x0, nv)
    nx, nv = x0,nv
    for k=1:100
        nx,nv = step(model, x0, nv)
    end
    return nx, nv
end


solve(model, x0, nv)

using Plots

kv = [e[1] for e in model.grid[1,:]]
plot(kv,[e[1] for e in nx[1,:]])
plot!(kv,[e[1] for e in nx[2,:]])
plot!(kv,kv*model.calibration.p.δ)
scatter!(model.calibration.s, model.calibration.x)

plot!(kv,[e[1] for e in nv[1,:]])
plot!(kv,[e[1] for e in nv[2,:]])


fobj = u->-Q(model, m, s, SVector(u...), φv)

using ForwardDiff

ForwardDiff.hessian(fobj, x)