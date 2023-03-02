using NoLib
using Plots

import NoLib: transition, arbitrage, recalibrate, initial_guess, projection, equilibrium

include("models/consumption_savings_ayiagari.jl")

# check we can solve the model with default calibration

@time sol = NoLib.time_iteration(model; verbose=true, improve=true, T=20)
x0 = sol.solution
μ0 = NoLib.ergodic_distribution(model, x0)



svec = [e[1] for e in model.grid[1,:]]

pl1 = plot(svec, [e[1] for e in x0[1,:]])
plot!(svec, [e[1] for e in x0[2,:]])
plot!(svec, [e[1] for e in x0[3,:]])

pl2 = plot(svec, [e[1] for e in μ0[1,:]])
plot!(svec, [e[1] for e in μ0[2,:]])
plot!(svec, [e[1] for e in μ0[3,:]])

plot(pl1, pl2)

### find the equilibrium

y0 = model.calibration.y

@time NoLib.NoLib.ss_residual(model, y0; x0=x0);


@time NoLib.NoLib.ss_residual(model, y0; x0=x0, diff=true);

@time ys = NoLib.find_equilibrium(model)



### visualize equilibrium

# plot residuals
kvec = range(30, 60; length=20)
rvec = []

@time NoLib.ss_residual(model, SVector(41.0); diff=true, x0=sol.solution)

using ProgressLogging
@progress for k in kvec
    val = try
        NoLib.ss_residual(model, SVector(k), improve=true; x0=sol.solution, diff=true)
    catch 
        NaN
    end
    push!(rvec, val)
end


#### check
## for some reason there are kinks in that graph

using Plots
plot(kvec, [e[1][1] for e in rvec])
for i=[ e for e in 1:length(rvec) if (e%5==0)]
    x_ = kvec[i]
    y_, c_ = rvec[i]
    scatter!([x_], [y_], color="red")
    plot!(kvec, y_ .+ (kvec.-x_).*c_, color="red", alpha=0.5)
end
scatter!(ys, [0], color="black")
# plot!(kvec, rvec2)
# scatter!(kvec, rvec)
plot!(kvec, kvec*0)


###
### Sequence-space jacobian

J = NoLib.f_residuals(model, μ0, x0, y0; diff=true);

## TEST Gauss elimination
# NoLib.gauss_elimination_residuals(R, J)
R = [μ0, x0, y0];

function build_X(y, J)
    T_x = J.F._x_1 \ J.F._x_2;

    # If y contains only 1 element
    if length(y) == 1
        interm = NoLib.invert(J.F._y * SVector(y), T_x)
        val = J.A._y * y + J.A._μ * NoLib.invert_G_μ(J.G._y * SVector(y), J.G._μ) + J.A._x * interm + J.A._μ * NoLib.invert_G_μ(J.G._x * interm, J.G._μ)
        return cat(val; dims=2)
    # If y has a higher dimension
    else
        X = Array{Any}(undef, length(y))
        for (i, e) in enumerate(y)
            interm = NoLib.invert(J.F._y * e, T_x);
            val = J.A._y * e + J.A._μ * NoLib.invert_G_μ(J.G._y * e, J.G._μ) + J.A._x * interm + J.A._μ * NoLib.invert_G_μ(J.G._x * interm, J.G._μ)    
            X[i] = val
        end
        X = cat(X; dims=2)
        return X
    end
end

rx = R[1];
ry = R[2];
rz = R[3];

T_x = J.F._x_1 \ J.F._x_2;

rx = NoLib.invert_G_μ(rx, J.G._μ);
rz = rz - J.A._μ * rx;
ry = NoLib.invert(ry, T_x);
rz = rz - (J.A._x * ry) + J.A._μ * NoLib.invert_G_μ(J.G._x * ry, J.G._μ);
X = build_X(rz, J);
rz = inv(X) * rz;
rz = SVector{1}(rz);
ry = ry + NoLib.invert(J.F._y * rz, T_x);
rx = rx + NoLib.invert_G_μ(J.G._y * rz, J.G._μ);
rx = rx + NoLib.invert_G_μ(J.G._x * ry, J.G._μ);


