module NoLib

    import Base: eltype
    using LabelledArrays
    using StaticArrays
    using ForwardDiff
    
    include("interp.jl")
    using NoLib.Interpolation: interp
    
    ⟂(a,b) = min(a,b)
    # ⟂ᶠ(a,b)

    import Base: getindex
    
    include("grids.jl")
    include("garray.jl")
    include("model.jl")
    include("simul.jl")
    include("time_iteration.jl")
    include("hetag.jl")

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