
has_complementarities(model) = false

function complementarities(model, s_::Tuple, x_::SVector)

    m,s = NoLib.split_states(model, s_)
    mm_ = NoLib.LVectorLike(model.calibration.m, m)
    ss_ = NoLib.LVectorLike(model.calibration.s, s)
    ss = catl(mm_, ss_)

    xx = NoLib.LVectorLike(model.calibration.x, x_)    

    complementarities(model, ss, xx)

end


complementarities(model, s::SLArray, x::SLArray) = error("Not Implemented")


function F(model, s, x::SVector, φ::Union{GArray, DFun})
    r = sum(
         w*arbitrage(model,s,x,S,φ(S)) 
         for (w,S) in τ(model, s, x)
    )
    if has_complementarities(model)
        comp = complementarities(model, s, x)
        r = r ⟂ comp
    end
    r
end


F(model, controls::GArray, φ::Union{GArray, DFun}) =
    GArray(
        model.grid,
        [
            F(model,s,x,φ)
            for (s,x) in zip(enum(model.grid), controls)
        ],
    )

# using LoopVectorization
function F!(out, model, controls::GArray, φ::Union{GArray, DFun})
    for n in 1:length(model.grid)
        i, j = NoLib.from_linear(model.grid, n)

        s_ = model.grid[n]
        s = ((i,j), s_)
        x = controls[n]
        out[n] = F(model,s,x,φ)
    end
end



# function F!(out, model, controls::GArray, φ::GArray)
#     @floop for (n, (s,x)) in enumerate(zip(enum(model.grid), controls))
#         out[n] = F(model,s,x,φ)
#     end
# end   #### no alloc

## no alloc
dF_1(model, s, x, φ) = ForwardDiff.jacobian(u->F(model, s, u, φ), x)

dF_1(model, controls::GArray, φ::Union{GArray, DFun}) =
    GArray(    # this shouldn't be needed
        model.grid,
        [
            dF_1(model, s, x, φ)
            for (s,x) in zip(enum(model.grid), controls) 
        ]
    )

function dF_1!(out, model, controls::GArray, φ::Union{GArray, DFun})
    for (n, (s,x)) in enumerate(zip(enum(model.grid), controls))
        out[n] = ForwardDiff.jacobian(u->F(model, s, u, φ), x)
    end
end    #### no alloc
    




dF_2(model, s, x::SVector, φ::GArray, dφ::GArray) = 
    sum(
            w*ForwardDiff.jacobian(u->arbitrage(model,s,x,S,u), φ(S))* dφ(S)
            for (w, S) in τ(model, s, x)
    )   ### no alloc


dF_2(model, controls::GArray, φ::GArray, dφ::GArray) =
    GArray(
        model.grid,
        [(dF_2(model,s,x,φ,dφ)) for (s,x) in zip(enum(model.grid), controls) ],
    )

function dF_2!(out::GArray, model, controls::GArray, φ::GArray, dφ::GArray)
    for (n, (s,x)) in enumerate(zip(enum(model.grid), controls))
        out[n] = dF_2(model, s, x, φ, dφ)
    end
end   

include("dev_L2.jl")

using LinearMaps


function time_iteration_workspace(model; interp_mode=:linear)

    x0 = (NoLib.initial_guess(model))
    x1 = deepcopy(x0)
    x2 = deepcopy(x0)
    r0 = deepcopy(x0)
    dx = deepcopy(x0)
    N = length(dx)
    n = length(dx.data[1])
    J = GArray(
        model.grid,
        zeros(SMatrix{n,n,Float64,n*n}, N)
    )
    φ = DFun(model, x0; interp_mode=interp_mode)
    return (;x0, x1, x2, r0, dx, J, φ)
end

function newton_workspace(model; interp_mode=:linear)

    
    res = time_iteration_workspace(model; interp_mode=interp_mode)
    T =  NoLib.dF_2(model, res.x0, res.φ)
    res = merge(res, (;T=T,memn=(;du=deepcopy(res.x0), dv=deepcopy(res.x0))))
    return res
end


function time_iteration(model, workspace=time_iteration_workspace(model);
    T=500, K=10, tol_ε=1e-8, tol_η=1e-6, verbose=false, improve=false, interp_mode=:cubic, engine=:none
    )

    # mem = typeof(workspace) <: Nothing ? time_iteration_workspace(model) : workspace
    mbsteps = 5
    lam = 0.5

    (;x0, x1, x2, dx, r0, J, φ) = workspace

    for t=1:T
        
        NoLib.fit!(φ, x0)

        if engine==:cpu
            F!(r0, model, x0, φ, CPU())
        else
            F!(r0, model, x0, φ)
        end
        # r0 = F(model, x0, φ)

        ε = norm(r0)

        if ε<tol_ε
            return (;message="Solution found", solution=x0, n_iterations=t, dr=φ)
        end

        # solve u->F(u,x0) 
        # result in x1
        # J and r0 are modified

        x1.data .= x0.data

        for k=1:10

            if engine==:cpu
                F!(r0, model, x1,  φ, CPU())
            else
                F!(r0, model, x1,  φ)
            end

            ε_n = norm(r0)
            if ε_n<tol_ε
                break
            end

            if engine==:cpu
                dF_1!(J, model, x1,  φ, CPU())
            else
                dF_1!(J, model, x1,  φ)
            end

            dx.data .= J.data .\ r0.data

            for k=0:mbsteps
                x2.data .= x1.data .- dx.data .* lam^k
                if engine==:cpu
                    F!(r0, model, x2,  φ, CPU())
                else
                    F!(r0, model, x2,  φ)
                end
                ε_b = norm(r0)
                if ε_b<ε_n
                    break
                end
            end

            x1.data .= x2.data

        end

        η = distance(x0, x1)
        verbose ? println("$t: $ε : $η: ") : nothing

        if !(improve)

            x0.data .= x1.data
            
        else
            # x = T(x)
            # x1 = T(x0)
            # x - x1 = -T'(x) (x - x0)
            # x = x1 - T' (x - x0)
            # x = (I-T')\(x1 - T' x0)

            J_1 = NoLib.dF_1(model, x1, φ)
            J_2 =  NoLib.dF_2(model, x1, φ)
            J_2.M_ij[:] *= -1.0
            Tp = J_1 \ J_2
            r = x1 - Tp * x0

            x0 = neumann(Tp, r; K=1000)
            x1.data .= x0.data

        end


    end

    return (;solution=x0, message="No Convergence") # The only allocation when workspace is preallocated

end


function newton(model, workspace=newton_workspace(model);
    K=10, tol_ε=1e-8, tol_η=1e-6, verbose=false, improve=false, interp_mode=:cubic
    )

    # mem = typeof(workspace) <: Nothing ? time_iteration_workspace(model) : workspace

    (;x0, x1, x2, dx, r0, J, φ, T, memn) = workspace


    for t=1:K
        
        NoLib.fit!(φ, x0)

        F!(r0, model, x0, φ)

        ε = norm(r0)
        verbose ? println("$t: $ε") : nothing

        if ε<tol_ε
            return (;message="Solution found", solution=x0, n_iterations=t, dr=φ)
        end

        x1.data .= x0.data


        dF_1!(J, model, x0, φ)
        dF_2!(T, model, x0, φ)

        
        T.M_ij .*= -1.0
        T.M_ij .= J.data .\ T.M_ij 
        
        r0.data .= J.data .\ r0.data
        
        neumann!(dx, T, r0, memn; K=1000)

        x0.data .= x1.data .- dx.data
        x1.data .= x0.data


        # end


    end

    return (;solution=x0, message="No Convergence") # The only allocation when workspace is preallocated

end