import Base: *

struct LF{G, T, n_s}
    grid::G
    M_ij::Matrix{T}
    S_ij::Matrix{Tuple{Tuple{Int64}, SVector{n_s, Float64}}}
end

function dF_2(model, x::GArray, φ::GArray )
    res = []
    for (s,x) in zip(enum(model.grid), x)
        l = []
        for (w, S) in τ(model, s, x)
            el = w*ForwardDiff.jacobian(u->arbitrage(model,s,x,S,u), φ(S))
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


function apply_L_2!(dr, L2, dφ)
    (;M_ij, S_ij) = L2
    N,K = size(M_ij)
    for n=1:N
        t0 = dr[n]*0.0
        for k=1:K
            F_x = M_ij[n,k]
            S = S_ij[n,k]
            t0 += F_x*dφ(S)
        end
        dr[n] = t0
    end
end

function apply_L_2(L2, dφ)
    dr = deepcopy(dφ)
    apply_L_2!(dr, L2, dφ)
    return dr
end

# function invert!(dx, dr, L2, dφ; K=1000)
#     # modifies dx and dr
#     for k=1:K
#         dx.data .+= dr.data
#         apply_L_2!(dr, L2, dφ)
#     end
# end


function invert(dr, L2; K=1000, τ_η=1e-10)
    dx = deepcopy(dr)
    du = deepcopy(dr)
    dv = deepcopy(dr)
    for k=1:K
        apply_L_2!(du, L2, dv)
        dx.data .+= du.data
        η = norm(du)
        if η<τ_η
            break
        end
        du,dv=dv,du
    end
    return dx
end

import Base: /,*
import LinearAlgebra: mul!

*(L::LF, x0) = apply_L_2( L, x0)
mul!(y0, L::LF, x0) = apply_L_2!(y0, L, x0)

\(J, L::LF) = LF(L.grid, J.data .\ L.M_ij, L.S_ij)
\(L::LF, r) = invert(r, L)

using LinearMaps
using BlockDiagonals
convert(::Type{BlockDiagonal}, ga::GArray{G, Vector{T}}) where T where G = BlockDiagonal(ga.data)


convert(::Type{Matrix}, L::LF) = convert(Matrix, convert(LinearMap,L))
convert(::Type{LinearMap}, L::LF) = let
    elt = eltype(L.M_ij)
    p,q = size(elt)
    N = length(L.grid)
    typ = SVector{q, Float64}
    fun = u->(ravel(L*GArray(L.grid, reinterpret(typ, u))))
    LinearMap(
        fun,
        p*N,
        q*N
    )
end



