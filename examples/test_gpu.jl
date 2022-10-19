using CUDA
using Test
using CUDA

# N = 100
# x_d = CUDA.fill(1.0f0, N)  # a vector stored on the GPU filled with 1.0 (Float32)
# y_d = CUDA.fill(2.0f0, N)  # a vector stored on the GPU filled with 2.0

# function gpu_add1!(y, x)
#     for i = 1:length(y)
#         @inbounds y[i] += x[i]
#     end
#     return nothing
# end

# fill!(y_d, 2)

# @cuda gpu_add1!(y_d, x_d)

# @test all(Array(y_d) .== 3.0f0)


# function gpu_add2!(y, x)
#     index = threadIdx().x    # this example only requires linear indexing, so just use `x`
#     stride = blockDim().x
#     for i = index:stride:length(y)
#         @inbounds y[i] += x[i]
#     end
#     return nothing
# end

# fill!(y_d, 2)
# @cuda threads=256 gpu_add2!(y_d, x_d)
# @test all(Array(y_d) .== 3.0f0)


# function gpu_add3!(model, y, x)
#     index = threadIdx().x    # this example only requires linear indexing, so just use `x`
#     stride = blockDim().x
#     st = sum(i for i=1:10)
    
#     for i = index:stride:length(y)
#         @inbounds y[i] += x[i] * st[1]
#     end
#     return nothing
# end


# fill!(y_d, 2)
# @cuda threads=256 gpu_add3!(model, y_d, x_d)
# @test all(Array(y_d) .== 3.0f0)


include("neoclassical_model.jl")

using CUDA


φ = GVector(model.grid, [Iterators.repeated(SVector(model.x), length(model.grid))...])
i0 = 3
s0_ = [NoLib.enum(model.grid)...][i0]
s0 = [NoLib.enum(model.grid)...][i0]
m0 = s0[2:end]

s0 = [enum(model.grid)...][1]
x0 = φ[1]

# function gpppu_!(model, y, r, φ)
#     index = threadIdx().x    # this example only requires linear indexing, so just use `x`
#     stride = blockDim().x
    
#     for j = index:stride:length(y)
#         # @inbounds s_d = model.grid.g1[1]
#         # @inbounds s_c = model.grid.g2[j]
#         # s = SVector(s_d..., s_c...)
#         @inbounds s_ = model.grid[1,j]
#         s = ((1,j), s_)
#         x = φ.data[j]
#         rr = sum(
#             w*NoLib.arbitrage(model,s,x,S,φ(S)) 
#             for (w,S) in NoLib.τ(model, s, x)
#         )
        
#         r.data[j] = rr
#     end
#     return nothing
# end


import Adapt

import Adapt
Adapt.@adapt_structure NoLib.GArray

r_gpu = Adapt.adapt(CuArray, deepcopy(φ*0) )
φ_gpu = Adapt.adapt(CuArray, φ) 


# N = 100
# x = CUDA.fill(1.0f0, N)  # a vector stored on the GPU filled with 1.0 (Float32)
# y = CUDA.fill(2.0f0, N)  # a vector stored on the GPU filled with 2.0

# @cuda threads=64 gpppu_!(model, y, r_gpu, φ_gpu)




using CUDA
using CUDAKernels # Required to access CUDADevice
using KernelAbstractions


# @kernel function mul2_kernel(A)
#     I = @index(Global)
#     A[I] = 2 * A[I]
#   end

# A = CuArray(ones(1024, 1024))
# ev = mul2_kernel(CUDADevice(), 16)(A, ndrange=size(A))
# wait(ev)
# all(A .== 2.0)
import KernelAbstractions.Extras.LoopInfo: @unroll

@kernel function F_gpu(@Const(model), r, @Const(φ))

    c = @index(Global, Cartesian)
    i,j = c.I

    n = @index(Global)

    @inbounds s_ = model.grid[i,j]
    s = ((i,j), s_)
    @inbounds x = φ[i,j]
    
    rr = x*0
    for (w,S) in NoLib.τ(model, s, x)
        rr += w*NoLib.arbitrage(model,s,x,S,φ(S))
    end
    # rr = sum(
    #     w*NoLib.arbitrage(model,s,x,S,φ(S)) 
    #     for (w,S) in NoLib.τ(model, s, x)
    # )
        
    # r[i,j] = rr
    @inbounds r.data[i,j] = rr

        
end

ev = KGPU(CUDADevice(), 128)

event = ev(model, r_gpu, φ_gpu,  ndrange=(2,500))
wait(event)

# r_sol = Adapt.adapt(Array, r_gpu)
# φ_sol = Adapt.adapt(Array, φ_gpu)

# r_ref = NoLib.F(model, φ, φ)

# δ = r_ref - r_sol


function timing(K, φ0)
    r = deepcopy(φ0)
    for k=1:K
        NoLib.F!(r, model, φ0, φ0)
    end
    return r
end

res_cpu = timing(10000, φ)
@time timing(10000, φ)


function timing_gpu(K, φ0)
    r_gpu = Adapt.adapt(CuArray, φ0)
    φ_gpu = Adapt.adapt(CuArray, φ0)
    local event
    for k=1:K
        event = ev(model, r_gpu, φ_gpu,  ndrange=(2,500))
            # @synchronize
        # wait(event)
    end
    wait(event)
    r_sol = Adapt.adapt(Array, r_gpu)
    return r_sol
end

res_gpu = timing_gpu(10000, φ)

@time timing_gpu(10000, φ)

err = maximum(u->maximum(abs.(u)), r_cpu - r_gpu)