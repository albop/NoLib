using Test
using CUDA
using Adapt
using CUDAKernels # Required to access CUDADevice
using KernelAbstractions


include("neoclassical_model.jl")
Adapt.@adapt_structure NoLib.GArray


function timing(model; K=1000, method=:cpu, M=nothing)

    φ0 = GVector(model.grid, [Iterators.repeated(SVector(model.calibration.x), length(model.grid))...])
    r0 = deepcopy(φ0)
    r = deepcopy(φ0)
    local event
    if method==:gpu
        r_gpu = Adapt.adapt(CuArray, φ0)
        φ_gpu = Adapt.adapt(CuArray, φ0)
    end
    for k=1:K
        if method==:cpu
            event = fun_cpu(model, r, φ0,  ndrange=(500,2))
            # r0.data .+= r.data*0.0001
        elseif method==:hand
            event = fun_hand(model, r, φ0; M=M)
                # r0.data .+= r.data*0.0001
        elseif method==:normal
            NoLib.F!(r, model, φ0, φ0)
            # r0.data .+= r.data*0.0001
        elseif method==:gpu
            event = fun_gpu(model, r_gpu, φ_gpu,  ndrange=(500,2))
            # wait(event)

        end
        # r *= 0.0000001
        # φ0 += r
    end
    if method in (:cpu, :gpu)
        wait(event)
    end
    if method == :gpu
        r = Adapt.adapt(Array, r_gpu)
    end
    return r

end




include("temp.jl")
eval(code)
import KernelAbstractions as KA

fun_cpu = fun_(CPU(), 8)

φ0 = GVector(model.grid, [Iterators.repeated(SVector(model.calibration.x), length(model.grid))...])

r = deepcopy(φ0)

fun_cpu(model, r, φ0,  ndrange=(500,2))


