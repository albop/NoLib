using FiniteDiff

function τ(model, ss::Tuple, a::SVector, p0, p1)


    p = model.calibration.p
    P = model.transition
    Q = model.grid.g1.points

    ss = cover(p0, ss)

    (i,_),(s_) = ss # get current state values

    # TODO: replace following block by one nonallocating function
    k  = length(model.calibration.m)
    l = length(model.calibration.s)
    m = SVector((s_[i] for i=1:k)...)
    s = SVector((s_[i] for i=k+1:(k+l))...)

    it = (
        (
            P[i,j],
            (
                (j,),
                let
                    q = cover(p1,Q[j])
                    SVector(
                        q...,
                        transition(model, m, s, a, q, p)...
                    )
                end
            )
        )
        for j in 1:size(P, 2)
    )

    it

end

@resumable function τ_fit(model, ss::Tuple, a::SVector, p0; linear_index=false)

    p = model.calibration.p
    P = model.transition
    Q = model.grid.g1.points


    i = ss[1][1]


    n_m = length(model.calibration.m)

    ss = cover(p0, ss)

    (i,_),(s_) = ss # get current state values

    # TODO: replace following block by one nonallocating function
    k  = length(model.calibration.m)
    l = length(model.calibration.s)
    m = SVector((s_[i] for i=1:k)...)
    s = SVector((s_[i] for i=k+1:(k+l))...)

    for j in 1:size(P, 2)

        S = transition(model, m, s, a, Q[j], p)

        for (w, i_S) in trembling__hand(model.grid.g2, S)

            res = (
                P[i,j]*w,

                (
                    (linear_index ? to__linear_index(model.grid, (j,i_S)) : (j,i_S)),

                    SVector(Q[j]..., model.grid.g2[i_S]...)
                )
            )
            if linear_index
                res
                # res::Tuple{Float64, Tuple{Int64, Tuple{SVector{1},SVector{1}}}}
            else
                res
                # res::Tuple{Float64, Tuple{Tuple{Int64, Int64}, Tuple{SVector{1},SVector{1}}}}
            end
            @yield res
        end
    end
end


function G(model, μ::GDist{T}, x, p0; diff=true) where T

    μ1 = GArray(μ.grid, zeros(Float64, length(μ)))
    for ss in enum(model.grid)
        ss = cover(p0,ss)
        a = x[ss[1]...]
        for (w, (ind, _)) in τ_fit(model, ss, a, p0)
            μ1[ind...] += w*μ[ind...]
        end
    end

    if !diff
        return μ1
    else
        P = transition_matrix(model, x)
        G_x = FiniteDiff.finite_difference_jacobian(
            u->ravel(G(model, μ, unravel(x,u), p0; diff=false)),
            ravel(x)
        )
        G_p = FiniteDiff.finite_difference_jacobian(
            u->ravel(G(model, μ, x, u; diff=false)),
            p0
        )
        
        return μ1, P, G_x, G_p
    end

end


cover(p, s::Tuple{IND, Tuple{T_m, T_s}})  where IND where T_m <: SVector{U,T} where T_s where U where T = (s[1], (cover(p, s[2][1]), s[2][2]))
cover(p, s::Tuple{IND, Tuple{TT}})  where IND where TT <: SVector{U,T} where U where T = (s[1], (cover(p, s[2])))
cover(p, s::Tuple{IND, TT})  where IND where TT <: SVector{U,T} where U where T = (s[1], (cover(p, s[2])))
# cover(p, s::Tuple{IND, TT})  where IND where TT where d = (s[1], (cover(p, s[2])))


F(model, s, x::SVector, φ::GArray, p0, p1) =
    sum(
        w*arbitrage(model,cover(p0,s),x,cover(p1,S),φ(S)) 
         for (w,S) in τ(model, s, x, p0, p1)
    )


function F(model, controls::GArray, φ::GArray, p0, p1; diff=false)

    r = GArray(
        model.grid,
        [
            F(model,s,x,φ,p0,p1)
            for (s,x) in zip(enum(model.grid), controls)
        ],
    )

    if !diff
        return r
    end
    
    r_1 = dF_1(model, controls, φ, p0, p1)
    r_2 = dF_2(model, controls, φ, p0, p1)

    N = length(model.grid)
    n_x = length(eltype(controls))
    n_p = length(p0)
    
    # ordering of the jacobian is: (n_x*N)*n_p
    r_p0 = FiniteDiff.finite_difference_jacobian(
        u->ravel(F(model, controls, φ, u, p1; diff=false)),
        p0
    )
    r_p1 = FiniteDiff.finite_difference_jacobian(
        u->ravel(F(model, controls, φ, p0, u; diff=false)),
        p1
    )

    rm_p0 = reshape( view(r_p0, :), (n_x, N, n_p) )
    rm_p1 = reshape( view(r_p1, :), (n_x, N, n_p) )

    rr_p0 = permutedims(rm_p0, (1,3,2))
    rr_p1 = permutedims(rm_p1, (1,3,2))

    elt = SMatrix{n_x, n_p, Float64, n_x*n_p}
    r_3 = GArray(model.grid, reinterpret(elt, rr_p0[:]))
    r_4 = GArray(model.grid, reinterpret(elt, rr_p1[:]))

    return r, r_1, r_2, r_3, r_4

end

dF_1(model, controls::GArray, φ::GArray, p0, p1) =
    GArray(    # this shouldn't be needed
        model.grid,
        [
            ForwardDiff.jacobian(u->F(model, s, u, φ, p0, p1), x)
            for (s,x) in zip(enum(model.grid), controls) 
        ]
    )

dF_2(model, s, x::SVector, φ::GArray, dφ::GArray, p0, p1) = 
    sum(
            w*ForwardDiff.jacobian(u->arbitrage(model,cover(p0,s),x,S,u), φ(S))* dφ(S)
            for (w, S) in τ(model, s, x, p0, p1)
    )


dF_2(model, controls::GArray, φ::GArray, dφ::GArray, p0, p1) =
    GArray(
        model.grid,
        [(dF_2(model,s,x,φ,dφ,p0,p1)) for (s,x) in zip(enum(model.grid), controls) ],
    )

using LinearMaps

# dF_2(model, x::GArray, φ::GArray, p0, p1) = let
#     n = length(φ)*length(eltype(φ))
#     fun = u->(ravel(dF_2(model, x, φ, unravel(x, u), p0, p1)))
#     xxx = ravel(x)
#     # return xxx
#     e = fun(xxx)
#     LinearMap(
#         fun,
#         n,
#         n
#     )
# end

function dF_2(model, x::GArray, φ::GArray, p0, p1 )
    # TODO: one should preallocate here
    res = []
    for (s,x) in zip(enum(model.grid), x)
        l = []
        for (w, S) in τ(model, cover(p0, s), x)
            el = -w*ForwardDiff.jacobian(u->arbitrage(model,cover(p0,s),x,cover(p1,S),u), φ(S))
            push!(l, (el, S))
        end
        push!(res, l)
    end
    N = length(res)
    J = length(res[1])
    tt = res[1][1]
    M_ij = Array{typeof(tt[1])}(undef,N,J)
    S_ij = Array{typeof(tt[2])}(undef,N,J)
    for n=1:N
        for j=1:J
            M_ij[n,j] = res[n][j][1]
            S_ij[n,j] = res[n][j][2]
        end
    end
    return LF(model.grid, M_ij, S_ij)
end
