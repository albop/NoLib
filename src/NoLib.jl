module NoLib

    import Base: eltype
    using LabelledArrays
    using StaticArrays
    using ForwardDiff


    include("interp.jl")
    
    using NoLib.Interpolation: interp
    import Base: getindex
    
    ⟂(a,b) = min(a,b)
    # ⟂ᶠ(a,b)


    abstract type AGrid{d} end

    struct CGrid{d} <: AGrid{d}
        ranges::NTuple{d, Tuple{Float64, Float64, Int64}}
    end

    getindex(g::CGrid{1}, i::Int) = SVector{1}(
        g.ranges[1][1] + (g.ranges[1][2]-g.ranges[1][1])*( (i-1)/(g.ranges[1][3]-1))
    )

    getindex(g::CGrid{1}, ::Colon) = [SVector(i) for i in range(g.ranges[1]...)]


    struct SGrid{d} <: AGrid{d}
        points::Vector{SVector{d, Float64}}
    end

    function SGrid(Q::Matrix)
        d = size(Q,2)
        return SGrid{d}([SVector(Q[i,:]...) for i=1:size(Q,1)])
    end

    struct PGrid{G1, G2, d} <: AGrid{d}
        g1::G1
        g2::G2
        points::Vector{SVector{d, Float64}}
    end

    getindex(g::SGrid{d}, ::Colon) where d = g.points
    getindex(g::SGrid{d}, i::Int) where d = g.points[i]

    cover(m,v::SVector{d,T}) where d where T = SVector{d,T}(m...,v[length(m)+1:end]...)

    PGrid(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid{SGrid{d1}, CGrid{d2}, d1+d2}(g1, g2,
        [SVector(v1...,v2...) for ((i,v1), (j,v2)) in Base.Iterators.product( iti( g1), iti(g2) ) ][:]
    )

    ×(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid(g1,g2)

    import Base: getindex
    getindex(g::PGrid{G1, G2, d}, i) where G1 where G2 where d = g.points[i]
    getindex(g::PGrid{G1, G2, d}, i::Int64, j::Int64) where G1 where G2 where d = g.points[ i + length(g.g1)*(j-1)]
    getindex(g::PGrid{G1, G2, d}, i::Int64, ::Colon) where G1 where G2 where d = g.g2[:] # TODO: should error if i out of bounds
    getindex(g::PGrid{G1, G2, d}, ::Colon, i::Int64) where G1 where G2 where d = g.g1[:]


    @inline to__linear_index(g::PGrid, ind::Tuple{Int64, Int64}) =  ind[1] + length(g.g1)*(ind[2]-1)


    import Base: iterate
    import Base: length
    import Base: getindex
    import Base: setindex!

    iti(sg::SGrid{d}) where d = enumerate(p for p in sg.points)
    iti(cg::CGrid{d}) where d = enumerate(  SVector(el...)    for el in Base.Iterators.product((range(r[1], r[2], r[3]) for r in cg.ranges)...) ) 
    iti(pg::PGrid) = ( ((i,j),(v1,v2)) for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )

    function enum(pg::PGrid; linear_index=false) 
        if !linear_index
            ( ( (i,j), SVector(v1...,v2...) )
              for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )

        else
            ( ( to__linear_index(pg,(i,j)), SVector(v1...,v2...) )
              for ((i,v1), (j,v2)) in Base.Iterators.product( (iti( pg.g1)), (iti(pg.g2)) ) )        end
    end

    length(pg::PGrid{G1, G2, d}) where G1 where G2 where d = length(pg.g1)*length(pg.g2)
    length(sg::SGrid{d}) where d = length(sg.points)
    length(cg::CGrid{d}) where d = prod(e[3] for e in cg.ranges)

    ## GArrays

    struct GArray{G,T}
        grid::G
        data::Vector{T}
    end

    GArray(grid::G, x::AbstractVector{T}) where G where T = GArray{G,T}(grid, copy(x))


    norm(v::GArray) = maximum(u->maximum(abs, u), v.data)

    import Base: size
    size(a::GArray{PGrid{G1, G2, d}, T}) where G1 where G2 where d where T = (length(a.grid.g1), length(a.grid.g2))


    # TODO: check
    getindex(a::GArray{PGrid{G1, G2, d}, T}, i::Int, j::Int) where G1 where G2 where d where T = a.data[ i + length(a.grid.g1)*(j-1)]
    getindex(a::GArray{PGrid{G1, G2, d}, T}, i::Int, ::Colon) where G1 where G2 where d where T = [a[i,j] for j=1:length(a.grid.g2)]
    getindex(a::GArray{PGrid{G1, G2, d}, T}, ::Colon, j::Int) where G1 where G2 where d where T = [a[i,j] for i=1:length(a.grid.g1)]

    getindex(a::GArray{PGrid{G1, G2, d}, T}, ::Colon) where G1 where G2 where d where T = a.data


    # TODO: check
    setindex!(a::GArray{PGrid{G1, G2, d}, T}, v, i::Int, j::Int) where G1 where G2 where d where T = (a.data[ i + length(a.grid.g1)*(j-1)] = v)



    GDist{G} = GArray{G, Float64}

    eltype(g::GArray{G,T}) where G where T = T

    # warning: these functions don't copy any data
    ravel(g::GArray) = reinterpret(Float64, g.data)
    unravel(g::GArray, x) = GArray(
        g.grid,
        reinterpret(eltype(g), x)
    )

    iterate(g::GArray) = iterate(g.data)
    iterate(g::GArray, i) = iterate(g.data, i)
    
    length(g::GArray) = length(g.data)
    getindex(g::GArray{G,T}, i::Int64) where G where T = g.data[i]
    setindex!(g::GArray, x, i) = (g.data[i] = x)
    

    # interpolating indexing
    function (xa::GArray{PGrid{G1, G2, d}, T})(i::Int64, p::SVector{d2, U}) where G1<:SGrid where G2<:CGrid{d2} where d where d2 where T where U
        g1 = xa.grid.g1
        g2 = xa.grid.g2
        dims = tuple(length(g1), (e[3] for e in g2.ranges)... )
        # ranges = tuple( (range(e...) for e in g2.ranges)... )
        v = view(reshape(xa.data, dims),i,:)
        res = interp(g2.ranges, v, p...)
        res
    end

    (xa::GArray{PGrid{G1, G2, d}, T})(i::Int64, j::Int64) where G1 where G2 where U where d where T  = xa[i,j]
    (xa::GArray{PGrid{G1, G2, d}, T})(S::Tuple{Tuple{Int64}, U}) where G1 where G2 where U where d where T = xa(S[1][1],S[2][2])
    (xa::GArray{PGrid{G1, G2, d}, T})(S::Tuple{Tuple{Int64, Int64}, U}) where G1 where G2 where U where d where T = xa[S[1]...]




    # enum(g::GArray) = enumerate(g)

    import Base: *, \, +, -, /

    *(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.*B.data)
    \(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.\B.data)
    /(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data./B.data)
    +(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.+B.data)
    -(A::GArray{G,T}, B::GArray{G,U}) where G where T  where U = GArray(A.grid, A.data.-B.data)
    *(A::GArray{G,T}, x::Number) where G where T  where U = GArray(A.grid, A.data .* x)
    
    *(A::GArray{G,T}, B::AbstractArray{Float64}) where G where T <:SMatrix{p, q, Float64, n}  where p where q where n = 
        ravel(
            GArray(
                A.grid,
                A.data .* reinterpret(SVector{q, Float64}, B)
            )
        )

    import Base: convert

    function Base.convert(::typeof(Matrix), A::GArray{G,T}) where G where T <:SMatrix{p, q, Float64, k}  where p where q where k
        N = length(A.data)
        n0 = N*p
        n1 = N*q
        M = zeros(n0, n1)
        for n=1:N
            i = p*(n-1)+1
            j = q*(n-1)+1
            M[i:i+p-1, j:j+q-1] = A.data[n]
        end
        return M
    end


    ## model functions

    import Base: merge
    function merge(a::SLArray, b::SLArray)
        # TODO: type is incorrect (should correspond to SVector)
        m = (merge(convert(NamedTuple,a), convert(NamedTuple,b)))
        return SLVector(m)
    end

    function LVectorLike(m0::SLArray{Tuple{d}, T, 1, d, M}, m) where d where T where M
        tt = eltype(m)
        TT = SLArray{Tuple{d}, tt, 1, d, M}
        return TT(m...)
    end

    label_GArray(m, g::GArray) = GArray(g.grid, [LVectorLike(m, e) for e in g.data])

    function transition(model, m, s, x, M, p)
        m = LVectorLike(model.m,m)
        s = LVectorLike(model.s,s)
        x = LVectorLike(model.x,x)
        M = LVectorLike(model.m,M)
        S = transition(model,m,s,x,M,p)
        return SVector(S...)
    end
    
    function arbitrage(model, m, s, x, M, S, X, p)
        m = LVectorLike(model.m, m)  # this does not keep the original type
        s = LVectorLike(model.s, s)
        x = LVectorLike(model.x, x)
        M = LVectorLike(model.m, M)
        S = LVectorLike(model.s, S)
        X = LVectorLike(model.x, X)
        r = arbitrage(model,m,s,x,M,S,X,p)
        return SVector(r...)
    end

    function arbitrage(model, s, x, S, X)
        p = model.p
        arbitrage(model, s[2][1], s[2][2], x, S[2][1], S[2][2], X, p)
    end


    ### transition function


    function τ(model, ss::Tuple, a::SVector)


        p = model.p
        (i,_),(m, s) = ss # get current state values

        Q = model.grid.g1.points

        j = 1

        it = (
            (
                model.P[i,j],
                (
                    (j,),
                    (
                        Q[j],
                        transition(model, m, s, a, Q[j], p)
                    )
                )
            )
            for j in 1:size(model.P, 2)
        )

        it

    end


    function trembling__hand(g::CGrid{1}, xv)
        x = xv[1]
        r = g.ranges[1]
        u = (x-r[1])/(r[2]-r[1])
        n = r[3]

        i_ = floor(Int, u*(n-1))
        i_ = max(min(i_,n-2),0)
        λ = u*(n-1)-i_
    
        λ = min(max(λ, 0.0), 1.0)
        
        (
            (1-λ, i_+1),
            (λ, i_+2)
        )
    end
    
    using ResumableFunctions

    @resumable function τ_fit(model, ss::Tuple, a::SVector; linear_index=false)

        p = model.p

        i = ss[1][1]

        Q = model.grid.g1.points

        n_m = length(model.m)
        (m, s) = (ss[2][1:n_m], ss[2][n_m+1:end])

        for j in 1:size(model.P, 2)

            S = transition(model, m, s, a, Q[j], p)

            for (w, i_S) in trembling__hand(model.grid.g2, S)

                res = (
                    model.P[i,j]*w,

                    (
                        (linear_index ? to__linear_index(model.grid, (j,i_S)) : (j,i_S)),

                        (Q[j], model.grid.g2[i_S])
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

    function G(model, μ::GDist{T}, x) where T
        μ1 = GArray(μ.grid, zeros(Float64, length(μ)))
        for ss in iti(model.grid)
            # a = x(ss[1], ss[3])
            a = x[ss[1]...]
            for (w, (ind, _)) in τ_fit(model, ss, a)
                μ1[ind...] += w*μ[ind...]
            end
        end
        μ1
    end

    # TODO: write interpolation version of G


    function transition_matrix(model, x)

        N = length(x)
        P = zeros(N,N)

        for (ss,a) in zip(enum(model.grid; linear_index=false),x)
            i = to__linear_index(model.grid, ss[1])
            for (w, (j, _)) in τ_fit(model, ss, a; linear_index=true)
                P[i,j] = w
            end
        end
        P

    end


    using LinearAlgebra: I

    function ergodic_distribution(model, x::GArray)

        P = transition_matrix(model, x)

        PP = [ (P-I)[:,1:end-1] ;;  ones(size(P,1))]
        R = zeros(size(PP,1))
        R[end] = 1
        μ = PP'\R
        
        ergo = GArray(model.grid, μ)

        return ergo

    end


    ## TODO: some ways to plot the ergo dist...

    τ(model, ss::Tuple, φ) = τ(model, ss, φ(ss))

    ## TODO: some simulation

    # TODO
    # function simulate(model, ss::Tuple, φ; T=20)
    #     Typ = typeof(ss)
    #     sim = Typ[ss]
    #     for t = 1:T
    #     end
    # end

    include("time_iteration.jl")

end # module


# function (xa::NoLib.GArray{PGrid{G1, G2, d}, T})(p::SVector{U, d}) where G1 where G2 where d where T where U
#     error("Not Implemented")
# end


# function (xa::NoLib.GArray{PGrid{G1, G2, d}, T})(p::SVector{d, Float64}) where G1<:CGrid where G2<:CGrid where d where T
#     g1 = xa.grid.g1
#     g2 = xa.grid.g2
#     dims = tuple( (e[3] for e in g1.ranges)..., (e[3] for e in g2.ranges)...)
#     ranges = tuple( (range(e...) for e in g1.ranges)..., (range(e...) for e in g2.ranges)...)
#     v = reshape(xa.data, dims)


# end