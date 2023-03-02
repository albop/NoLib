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
GridSpace(v::Vector{SVector{d, Float64}}) where d where N = GridSpace{length(v), d, Val{(:i_)}}(SVector(v...))

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



using .splines
using .splines: Linear, Cubic, MLinear, MCubic, CubicInterpolator, SplineInterpolator



struct DFun{Dom, Gar, Itp, vars}
    domain::Dom
    values::Gar
    itp::Itp
end

vvars(::DFun{D,G,I,vars}) where D where G where I where vars<:Val{t} where t = t
vars(dfun::DFun) = vvars(dfun)

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
# works for cartesian only
function DFun(domain, values::GVector{G,V}; interp_mode=:linear) where V where G<:CGrid
    vv = reshape(values.data, (e[3] for e in values.grid.ranges)...)
    if interp_mode == :linear
        itp = SplineInterpolator(values.grid.ranges; values=vv, k=1)
    elseif interp_mode == :cubic
        itp = SplineInterpolator(values.grid.ranges;  values=vv,k=3)       
    else
        throw("Unkown interpolation mode $(interp_mode)")
    end
    return DFun(domain, values, itp)
end

function DFun(domain, values::GVector{G,V}; interp_mode=:linear) where V where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid
    if interp_mode == :linear
        k=1
    elseif interp_mode == :cubic
        k=3
    else
        throw("Unkown interpolation mode $(interp_mode)")
    end

    # TODO: check values.data[i,:]
    sz = (e[3] for e in values.grid.g2.ranges)
    itps = tuple( (SplineInterpolator(values.grid.g2.ranges;  values=reshape(values[i,:], sz...),k=3)  for i=1:length(values.grid.g1)  )...)
    return DFun(domain, values, itps)

end

function DFun(model::DModel, values::GVector{G,V}; interp_mode=:linear) where V where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid

    domain = model.domain
    
    if interp_mode == :linear
        k=1
    elseif interp_mode == :cubic
        k=3
    else
        throw("Unkown interpolation mode $(interp_mode)")
    end

    # TODO: check values.data[i,:]
    sz = (e[3] for e in values.grid.g2.ranges)
    itps = tuple( (SplineInterpolator(values.grid.g2.ranges;  values=reshape(values[i,:], sz...),k=3)  for i=1:length(values.grid.g1)  )...)
    return DFun(domain, values, itps)

end

## Cart
function (f::DFun{A,B,I,vars})(x::SVector{d2, U})  where A where B<:GArray{G,V} where V where I where G<:CGrid where vars where d2 where U
    f.itp(x)
end

## PGrid
function (f::DFun{A,B,I,vars})(i::Int, x::SVector{d2, U})  where A where B<:GArray{G,V} where V where I where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid where vars where d2 where U
    f.itp[i](x)
end

function (f::DFun{A,B,I,vars})(loc::Tuple{Tuple{Int64}, SVector{d2, U}})  where A where B<:GArray{G,V} where V where I where G<:PGrid{G1,G2} where G1<:SGrid where G2<:CGrid where vars where d2 where U
    # TODO: not beautiful
    x_ = loc[2]
    dd1 = ndims(f.values.grid.g1)
    dd2 = ndims(f.values.grid.g2)
    x = SVector((x_[i] for i=(dd1+1):(dd1+dd2))...)
    println(x)
    f.itp[loc[1][1]](x)
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
