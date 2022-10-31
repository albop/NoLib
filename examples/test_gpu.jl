using CUDA
using Test
using CUDA
using Adapt


include("neoclassical_model.jl")

Adapt.@adapt_structure NoLib.GArray

using CUDA


φ = GVector(model.grid, [Iterators.repeated(SVector(model.x), length(model.grid))...])
# i0 = 3
# s0_ = [NoLib.enum(model.grid)...][i0]
# s0 = [NoLib.enum(model.grid)...][i0]
# m0 = s0[2:end]

# s0 = [enum(model.grid)...][1]
# x0 = φ[1]

using CUDAKernels # Required to access CUDADevice
using KernelAbstractions

@kernel function F_(@Const(model), r, @Const(φ))

    c = @index(Global, Cartesian)
    j,i = c.I

    n = @index(Global)

    s_ = model.grid[i,j]
    s = ((i,j), s_)
    x = φ[i,j]
    
    rr = x*0
    for (w,S) in NoLib.τ(model, s, x)
        rr += w*NoLib.arbitrage(model,s,x,S,φ(S))
    end
    # rr = sum(
    #     w*NoLib.arbitrage(model,s,x,S,φ(S)) 
    #     for (w,S) in NoLib.τ(model, s, x)
    # )
        
    r[n] = rr

end

using FLoops

function F_hand(model, r, φ; M=M)

    N = length(r)

    @floop for c in CartesianIndices((500,2))

        j = c[1]
        i = c[2]

        s_ = model.grid[i,j]
        s = ((i,j), s_)
        x = φ[i,j]
    
        rr = x*0
        for (w,S) in NoLib.τ(model, s, x)
            rr += w*NoLib.arbitrage(model,s,x,S,φ(S))
        end

        r[i,j] = rr
    end

end

F_gpu = F_(CUDADevice(), 128)
F_cpu = F_(CPU(), 8)


function timing(model; K=1000, method=:cpu, M=nothing)

    φ0 = GVector(model.grid, [Iterators.repeated(SVector(model.x), length(model.grid))...])
    r0 = deepcopy(φ0)
    r = deepcopy(φ0)
    local event
    if method==:gpu
        r_gpu = Adapt.adapt(CuArray, φ0)
        φ_gpu = Adapt.adapt(CuArray, φ0)
    end
    for k=1:K
        if method==:cpu
            event = F_cpu(model, r, φ0,  ndrange=(500,2))
            # r0.data .+= r.data*0.0001
        elseif method==:hand
            event = F_hand(model, r, φ0; M=M)
                # r0.data .+= r.data*0.0001
        elseif method==:normal
            NoLib.F!(r, model, φ0, φ0)
            # r0.data .+= r.data*0.0001
        elseif method==:gpu
            event = F_gpu(model, r_gpu, φ_gpu,  ndrange=(500,2))
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

M = zeros(2,500)
res_normal = timing(model; K=1, method=:normal)

res_hand = timing(model; K=1, method=:hand, M=M)

res_cpu = timing(model; K=1)
res_gpu = timing(model; K=1, method=:gpu)


@time res_normal = timing(model; K=10000, method=:normal);

@time res_hand = timing(model; K=10000, method=:hand);

@time res_cpu = timing(model; K=10000);
@time res_gpu = timing(model; K=10000, method=:gpu);


err = maximum(u->maximum(abs.(u)), res_cpu - res_normal)
err = maximum(u->maximum(abs.(u)), res_gpu - res_normal)

err = maximum(u->maximum(abs.(u)), res_cpu - res_gpu)

rr = NoLib.F(model, φ, φ)

r = deepcopy(φ)
NoLib.F!(r, model, φ, φ)




########
########
quote 
function cpu_F_(__ctx__, model, r, φ; )
    let model = (KernelAbstractions.constify)(model), φ = (KernelAbstractions.constify)(φ)
        $(Expr(:aliasscope))
        begin
            var"##N#385" = length((KernelAbstractions.__workitems_iterspace)(__ctx__))
            begin
                #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:264 =#
                for var"##I#384" = (KernelAbstractions.__workitems_iterspace)(__ctx__)
                    #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:265 =#
                    (KernelAbstractions.__validindex)(__ctx__, var"##I#384") || continue
                    #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:266 =#
                    c = KernelAbstractions.__index_Global_Cartesian(__ctx__, var"##I#384")
                    n = KernelAbstractions.__index_Global_Linear(__ctx__, var"##I#384")
                    #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:267 =#
                    begin
                        #= REPL[27]:1 =#
                        #= REPL[27]:3 =#
                        #= REPL[27]:4 =#
                        (j, i) = c.I
                        #= REPL[27]:6 =#
                        #= REPL[27]:8 =#
                        s_ = model.grid[i, j]
                        #= REPL[27]:9 =#
                        s = ((i, j), s_)
                        #= REPL[27]:10 =#
                        x = φ[i, j]
                        #= REPL[27]:12 =#
                        rr = x * 0
                        #= REPL[27]:13 =#
                        for (w, S) = NoLib.τ(model, s, x)
                            #= REPL[27]:14 =#
                            rr += w * NoLib.arbitrage(model, s, x, S, φ(S))
                            #= REPL[27]:15 =#
                        end
                        #= REPL[27]:21 =#
                        r[n] = rr
                    end
                    #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:268 =#
                end
            end
        end
        $(Expr(:popaliasscope))
        return nothing
    end
end

begin
    #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:52 =#
    if !($(Expr(:isdefined, :F_)))
        #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:53 =#
        begin
            $(Expr(:meta, :doc))
            F_(dev) = begin
                    #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:53 =#
                    F_(dev, (KernelAbstractions.NDIteration.DynamicSize)(), (KernelAbstractions.NDIteration.DynamicSize)())
                end
        end
        #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:54 =#
        F_(dev, size) = begin
                #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:54 =#
                F_(dev, (KernelAbstractions.NDIteration.StaticSize)(size), (KernelAbstractions.NDIteration.DynamicSize)())
            end
        #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:55 =#
        F_(dev, size, range) = begin
                #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:55 =#
                F_(dev, (KernelAbstractions.NDIteration.StaticSize)(size), (KernelAbstractions.NDIteration.StaticSize)(range))
            end
        #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:56 =#
        function F_(dev::Dev, sz::S, range::NDRange) where {Dev, S <: KernelAbstractions.NDIteration._Size, NDRange <: KernelAbstractions.NDIteration._Size}
            #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:56 =#
            #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:57 =#
            if (KernelAbstractions.isgpu)(dev)
                #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:58 =#
                return (KernelAbstractions.construct)(dev, sz, range, gpu_F_)
            else
                #= /home/pablo/.julia/packages/KernelAbstractions/DqITC/src/macros.jl:60 =#
                return (KernelAbstractions.construct)(dev, sz, range, cpu_F_)
            end
        end
    end
end
end
