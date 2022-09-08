using KernelAbstractions
using StaticArrays
using NoLib: SGrid, CGrid, PGrid, GArray
include("rbc_def.jl")

s0 = [NoLib.iti(model.grid)...][1]
x0 = SVector(model.x)
xvec = [x0 for i=1:length(model.grid)]

xv = GArray(model.grid, xvec)

@time NoLib.F(model, s0, x0, xv)


using KernelAbstractions

import NoLib: F

@kernel function FK!(model, out, xv, xvt)
    I = @index(Global)
    ss = model.grid[I]
    s = (1, ss[1], ss[2])
    x = xv[I]
    r = F(model, s, x, xvt)
    out[I] = r
end


function doit(model, xv)
    out = xv*0
    kernel = FK!(CPU(), 1)
    event = kernel(model, out, xv, xv, ndrange=length(xv))
    wait(event)
    return out
end

using BenchmarkTools


@time doit(model, xv);
@time F(model, xv, xv);

@benchmark F(model, xv, xv)
@benchmark doit(model, xv)




function FK_test(i, model, out, xv, xvt)
    # I = @index(Global)
    # s = (1,model.grid[i])
    ss = model.grid[i]
    s = (1, ss[1], ss[2])
    x = xv[i]
    println(s, x)
    r = F(model, s, x, xvt)
    out[i] = r
end

FK_test(3, model, out, xv, xv)
