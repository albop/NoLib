using Test


include("neoclassical_model.jl")

using CUDA


φ = GVector(model.grid, [Iterators.repeated(SVector(model.calibration.x), length(model.grid))...])
i0 = 3
s0_ = [NoLib.enum(model.grid)...][i0]
s0 = [NoLib.enum(model.grid)...][i0]
m0 = s0[2:end]

s0 = [enum(model.grid)...][1]
x0 = φ[1]


import Adapt

Adapt.@adapt_structure NoLib.GArray

using KernelAbstractions

import KernelAbstractions.Extras.LoopInfo: @unroll

@kernel function F_gpu(@Const(model), r, @Const(φ))

    c = @index(Global, Cartesian)
    i,j = c.I

    n = @index(Global)

    s_ = model.grid[i,j]
    s = ((i,j), s_)
    x = φ[i,j]
    
    rr = x*0
    for (w,S) in NoLib.τ(model, s, x)
        rr += w*NoLib.arbitrage(model,s,x,S,φ(S))
    end
    r[i,j] = rr

        
end

ev = F_gpu(CUDADevice(), 128)

ev2 = F_gpu(CPU(), 4)


function timing(K, φ0)
    r = deepcopy(φ0)
    local event
    for k=1:K
        NoLib.F!(r, model, φ0, φ0)
        # event = ev2(model, r, φ0,  ndrange=(2,500))
        r *= 0.0000001
        φ0 += r
    end
    # wait(event)
    return r

end

res_cpu = timing(10000, φ)
@time timing(10000, φ);


function timing_gpu(K, φ0)
    r_gpu = Adapt.adapt(CuArray, φ0)
    φ_gpu = Adapt.adapt(CuArray, φ0)
    local event
    for k=1:K
        event = ev(model, r_gpu, φ_gpu,  ndrange=(2,500))
            # @synchronize
        # wait(event)
        r_gpu *= 0.0000001
        φ_gpu += r_gpu

    end
    wait(event)
    r_sol = Adapt.adapt(Array, r_gpu)
    return r_sol
end

res_gpu = timing_gpu(10000, φ)

@time timing_gpu(10000, φ);

err = maximum(u->maximum(abs.(u)), res_cpu - res_gpu)