import Base: dims

abstract type Space{d} end

struct CartesianSpace{d,dims}
    # names::NTuple{d, Symbol}
    min::NTuple{d, Float64}
    max::NTuple{d, Float64}
end

CartesianSpace(a::Tuple{Float64}, b::Tuple{Float64}) = CartesianSpace{length(a), Val{(:x,)}}(a,b)
CartesianSpace(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64}) = CartesianSpace{length(a), Val{(:x_1, :x_2)}}(a,b)

function CartesianSpace(;kwargs...)
    names = tuple(keys(kwargs)...)
    a = tuple((v[1] for v in values(kwargs))...)
    b = tuple((v[2] for v in values(kwargs))...)
    d = length(names)
    println(names)
    println(a,b)
    return CartesianSpace{d, Val{names}}(a,b)
end

# TODO: why is this not working?
# dims(dom::CartesianSpace{d,dims}) where d where dims = dims
ddims(dom::CartesianSpace{d,dims}) where d where dims<:Val{t} where t = t

dims(dom::CartesianSpace) = ddims(dom)

ndims(dom::CartesianSpace{d, dims}) where d where dims = d



struct GridSpace{N,d,dims}
    points::SVector{N,SVector{d,Float64}}
end
GridSpace(v::SVector{N, SVector{d, Float64}}) where d where N = GridSpace{length(v), d, Val{(:i_)}}(SVector(v...))

ndims(gd::GridSpace{N,d,dims}) where N where d where dims = d
ddims(gd::GridSpace{N,d,dims}) where N where d where dims<:Val{e} where e = begin println("JI"); e end
dims(gd::GridSpace) = ddims(gd)


struct ProductSpace{A,B}
    spaces::Tuple{A,B}
end

cross(A::DA, B::DB) where DA<:GridSpace where DB<:CartesianSpace = ProductSpace{DA, DB}((A,B))
ndims(p::P) where P<:ProductSpace = ndims(p.spaces[1]) + ndims(p.spaces[2])
dims(p::P) where P<:ProductSpace = tuple(dims(p.spaces[1])..., dims(p.spaces[2])...)



abstract type Function{d} end




struct DFun{Dom, Gar, Itp, vars}
    domain::Dom
    values::Gar
    itp::Itp
end

DFun(domain, values, itp, vars) = DFun{typeof(domain), typeof(values), typeof(itp), Val{vars}}(domain, values, itp)
function DFun(domain, values, itp)
    if eltype(values) <: Number
        vars = :y
    elseif eltype(values) <: SVector && length(eltype(values)) == 1
        vars = (:y,)
    end
    return DFun(domain, values, itp, vars)
end
# DFun(domain, values, itp) = begin @assert ndims(domain)==1 ; DFun(domain, values, itp, (:y,)) end
function DFun(domain, values; interp_mode=:linear)
    if interp_mode == :linear
        return DFun(domain, values, MLinear())
    elseif interp_mode == :cubic
        vv = reshape(values.data, (e[3] for e in values.grid.ranges)...)
        itp = CubicInterpolator(values.grid; values=vv)
        return DFun(domain, values, itp)
    else
        throw("Unkown interpolation mode $(interp_mode)")
    end
end

using .splines

# function DFun(domain, values, itp::Cubic, vars)

# end


vvars(::DFun{D,G,I,vars}) where D where G where I where vars<:Val{t} where t = t
vars(dfun::DFun) = vvars(dfun)
# DFun(domain<:Dom, values) where Dom<:Space{1} = DFun(domain, values, MLinear())

# Gar<:CartesianGrid
# 1d

function (f::DFun{A,B,I,vars})(i, x::SVector{d2, U})  where A where B<:GArray{G,V} where V where I<:Linear where G<:PGrid{g1, g2}  where g1<:SSGrid where g2<:CGrid where vars where d2 where U
    
    gg1 = f.values.grid.g1
    ranges = f.values.grid.g2.ranges
    data = f.values.data
    dims = tuple(length(gg1), (e[3] for e in ranges)... )
    v = reshape(view(data, :), dims) 
    return interp(ranges, view(v,i,:), x...)

end

function (f::DFun{A,B,I,vars})(x::SVector{d2, U})  where A where B<:GArray{G,V} where V where I<:Linear where G<:CGrid where vars where d2 where U
    ranges = f.values.grid.ranges
    data = f.values.data
    dims = tuple( (e[3] for e in ranges)... )
    v = reshape(view(data, :), dims) 
    interp(ranges, v, x...)
end

function (f::DFun{A,B,I,vars})(x::SVector{d2, U})  where A where B<:GArray{G,V} where V where I<:Cubic where G<:CGrid where vars where d2 where U

    a = [e[1] for e in f.values.grid.ranges]
    b = [e[2] for e in f.values.grid.ranges]
    n = [e[3] for e in f.values.grid.ranges]
    splines.eval_UC_spline(a,b,n, f.itp.θ, x)
end



# Compatibility calls

(f::DFun)(x::Float64) = f(SVector(x))
(f::DFun)(x::Float64, y::Float64) = f(SVector(x,y))
(f::DFun)(x::Vector{SVector{d,Float64}}) where d = [f(e) for e in x]


ndims(df::DFun) = ndims(df.domain)




##
## functions defined by hand
##

struct Fun{Dom, FF}
    domain::Dom
    f::FF
end

Fun(f) = Fun(
    begin n=methods(f)[1].nargs-1; CartesianSpace( tuple( (-Inf for i=1:n)...), tuple( (Inf for i=1:n)...) ) end,
    f
)

function (f::Fun{Dom})(args...) where Dom
    return f.f(args...)
end

# function ddFun(dom, gar)
#     d = ndims(dom)
#     dt = typeof(dom)
#     gt = typeof(gar)
#     DFun{d, dt, gt}(dom, gar)
# end

function discretize(space::D2; n=ntuple(u->15,d))  where  D2<:CartesianSpace{d} where d
    
    if false in (isfinite.(space.max) .& isfinite.(space.min))
        # TODO : improve
        throw(DomainError("Impossible to discretize space with infinite bounds."))
    end

    grid = CGrid(
        tuple( 
            ( (space.min[i], space.max[i], n[i]) for i=1:ndims(space) )...
        )
    )

    return grid

end


function discretize(dom::ProductSpace{D1,D2}; n=ntuple(u->15,ndims(dom))) where D1<:GridSpace where  D2<:CartesianSpace
    grid_a = SSGrid(dom.spaces[1].points)
    g_b = dom.spaces[2]
    grid_b = CGrid(
        tuple( 
            ( (g_b.min[i], g_b.max[i], n[i]) for i=1:ndims(g_b) )...
        )
    )
    return grid_a × grid_b
end
