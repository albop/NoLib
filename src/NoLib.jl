module NoLib

    import Base: eltype
    using LabelledArrays
    using StaticArrays
    using ForwardDiff


    include("interp.jl")
    
    using NoLib.Interpolation: interp
    import Base: getindex
    
    abstract type AGrid{d} end

    struct CGrid{d} <: AGrid{d}
        ranges::NTuple{d, Tuple{Float64, Float64, Int64}}
    end

    getindex(g::CGrid{1}, i::Int) = SVector{1}(
        g.ranges[1][1] + (g.ranges[1][2]-g.ranges[1][1])*( (i-1)/(g.ranges[1][3]-1))
    )

    struct SGrid{d} <: AGrid{d}
        points::Vector{SVector{d, Float64}}
    end

    struct PGrid{G1, G2, d} <: AGrid{d}
        g1::G1
        g2::G2
        points::Vector{SVector{d, Float64}}
    end

    PGrid(g1::SGrid{d1}, g2::CGrid{d2}) where d1 where d2 = PGrid{SGrid{d1}, CGrid{d2}, d1+d2}(g1, g2,
    [SVector(v1...,v2...) for (v1, v2) in Base.Iterators.product( iti( g1), iti(g2) ) ][:]
    )

    import Base: getindex
    getindex(g::PGrid{G1, G2, d}, i) where G1 where G2 where d = g.points[i]
    getindex(g::PGrid{G1, G2, d}, i::Int64, j::Int64) where G1 where G2 where d = g.points[ i + length(g.g1)*(j-1)]


    # TODO: check
    getindex(a::GArray{PGrid{G1, G2, d}, T}, i::Int, j::Int) where G1 where G2 where d where T = a.data[ i + length(a.grid.g1)*(j-1)]

    # TODO: check
    setindex!(a::GArray{PGrid{G1, G2, d}, T}, v, i::Int, j::Int) where G1 where G2 where d where T = (a.data[ i + length(a.grid.g1)*(j-1)] = v)


    import Base: iterate
    import Base: length
    import Base: getindex
    import Base: setindex!

    iti(sg::SGrid{d}) where d = (p for p in sg.points)
    iti(cg::CGrid{d}) where d = (  SVector(el...)    for el in Base.Iterators.product((range(r[1], r[2], r[3]) for r in cg.ranges)...) ) 
    iti(pg::PGrid) = ((i,v1,v2) for ((i,v1), v2) in Base.Iterators.product( enumerate(iti( pg.g1)), iti(pg.g2) ) )

    length(pg::PGrid{G1, G2, d}) where G1 where G2 where d = length(pg.g1)*length(pg.g2)
    length(sg::SGrid{d}) where d = length(sg.points)
    length(cg::CGrid{d}) where d = prod(e[3] for e in cg.ranges)


    function LVectorLike(m0::SLArray{Tuple{d}, T, 1, d, M}, m) where d where T where M
        tt = eltype(m)
        TT = SLArray{Tuple{d}, tt, 1, d, M}
        return TT(m...)
    end

    struct GArray{G,T}
        grid::G
        data::Vector{T}
    end

    GArray(grid::G, x::AbstractVector{T}) where G where T = GArray{G,T}(grid, copy(x))


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
        arbitrage(model, s[2], s[3], x, S[2], S[3], X, p)
    end

    function τ(model, ss::Tuple, a::SVector)


        p = model.p
        i, m, s = ss # get current state values

        it = (
            (
                model.P[i,j],
                (
                    j,
                    model.Q[j,:],
                    transition(model, m, s, a, model.Q[j,:], p)
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

    @resumable function τ_fit(model, ss::Tuple, a::SVector)

        p = model.p
        i, m, s = ss # get current state values

        for j in 1:size(model.P, 2)
            S = transition(model, m, s, a, model.Q[j,:], p)
            for (w, i_S) in trembling__hand(model.grid.g2, S)

                res = (
                    model.P[i,j]*w,
                    (
                        j,
                        model.Q[j,:],
                        i_S,
                        model.grid.g2[i_S]
                    )
                )
                res::Tuple{Float64, Tuple{Int64, SVector{1}, Int64, SVector{1}}}
                @yield res
            end
        end
    end


    τ(model, ss::Tuple, φ) = τ(model, ss, φ(ss[1],ss[3]))

    function G(model, μ::GDist{T}, x) where T
        μ1 = GArray(μ.grid, zeros(Float64, length(μ)))
        for ss in iti(model.grid)
            a = x(ss[1], ss[3])
            for (w,(i, m, j, s)) in τ_fit(model, ss, a)
                μ1[i,j] += w
            end
        end
        μ1
    end

    # TODO
    # function simulate(model, ss::Tuple, φ; T=20)
    #     Typ = typeof(ss)
    #     sim = Typ[ss]
    #     for t = 1:T
    #     end
    # end

    function F0(model, s, x::SVector, xfut::GArray)
        tot = SVector((x*0)...)
        for (w, S) in τ(model, s, x)
            ind = (S[1], S[3])
            X = xfut(ind...)
            tot += w*arbitrage(model,s,x,S,X)
        end
        return tot
    end


    F(model, s, x::SVector, φ::GArray) = 
        sum(
             w*arbitrage(model,s,x,S,φ(S[1], S[3])) 
             for (w, S) in τ(model, s, x)
        )

    # F(model, s, x::SVector, φ::GArray).dφ = 
    #     sum(
    #          w*arbitrage(model,s,x,S,φ(S)).dφ(S)
    #          for (w, S) in τ(model, s, x)
    #     )

    F(model, controls::GArray, φ::GArray) =
        GArray(
            model.grid,
            [F(model,s,x,φ) for (s,x) in zip(iti(model.grid), controls) ],
        )

    function F!(out, model, controls, φ) 
        # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
        n=0
        for s in iti(model.grid)
            n += 1
            x = controls.data[n]
            out.data[n] = F(model,s,x,φ)
        end
        # end
    end

    dF(model, controls::GArray, φ::GArray) =
        GArray(    # this shouldn't be needed
            model.grid,
            [
                ForwardDiff.jacobian(u->F(model, s, u, φ), x)
                for (s,x) in zip(iti(model.grid), controls) 
            ]
        )

    function dF!(out, model, controls, φ) 
        # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
        n=0
        for s in iti(model.grid)
            n += 1
            x = controls.data[n]
            out.data[n] = ForwardDiff.jacobian(u->F(model, s, u, φ), x)
        end
        # end
    end

    # function dF2!(out, model, controls, φ) 
    #     # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
    #     n=0
    #     for s in iti(model.grid)
    #         n += 1
    #         x = controls.data[n]
    #         out.data[n] = ForwardDiff.jacobian(u->F(model, s, x, φ), φ)
    #     end
    #     # end
    # end
    

    FdF(model, controls::GArray, φ::GArray) =
        GArray(
            model.grid,
            [
                (F(model,s,x,φ), ForwardDiff.jacobian(u->F(model, s, u, φ), x))
                for (s,x) in zip(iti(model.grid), controls) 
            ]
        )


    function F0(model, controls::GArray, xfut::GArray)

        N = length(controls)
        res = GArray(
            model.grid,
            zeros(typeof(controls[1]), N)
        )
        for (i,(s,x)) in enumerate(zip(iti(model.grid), controls))
            res[i] = F(model,s,x,xfut)
        end
        return res
    end



    function dF0(model, controls::GArray, xfut::GArray)

        N = length(controls)
        res = deepcopy(controls)
        dres = GArray(
            model.grid,
            zeros(typeof(res[1]*res[1]'), N)
        )
        for (i,(s,x)) in enumerate(zip(iti(model.grid), controls))
            res[i] = F(model,s,x,xfut)
            dres[i] = ForwardDiff.jacobian(u->F(model, s, u, xfut), x)
        end
        return res, dres
    end
    

    # interpolating indexing
    function (xa::GArray{PGrid{G1, G2, d}, T})(i::Int64, p::SVector{d2, U}) where G1<:SGrid where G2<:CGrid where d where d2 where T where U
        g1 = xa.grid.g1
        g2 = xa.grid.g2
        dims = tuple(length(g1), (e[3] for e in g2.ranges)... )
        # ranges = tuple( (range(e...) for e in g2.ranges)... )
        v = view(reshape(xa.data, dims),i,:)
        res = interp(g2.ranges, v, p...)
        res
    end
    

    dF2(model, s, x::SVector, φ::GArray, dφ::GArray) = 
        sum(
             w*ForwardDiff.jacobian(u->arbitrage(model,s,x,S,u), φ(S[1], S[3]))* dφ(S[1], S[3])
             for (w, S) in τ(model, s, x)
        )
    
    dF2(model, controls::GArray, φ::GArray, dφ::GArray) =
        GArray(
            model.grid,
            [dF2(model,s,x,φ,dφ) for (s,x) in zip(iti(model.grid), controls) ],
        )

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