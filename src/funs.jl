abstract type Space{d} end

struct CartesianSpace{d}
    # names::NTuple{d, Symbol}
    min::NTuple{d, Float64}
    max::NTuple{d, Float64}
end
ndims(cd::CartesianSpace{d}) where d = d

struct GridSpace{N,d}
    points::SVector{N,SVector{d,Float64}}
end
ndims(gd::GridSpace{d}) where d = d

GridSpace(v::Vector{SVector{d, Float64}}) where d = GridSpace(SVector(v...))

struct ProductSpace{d,A,B}
    grid_A::A
    grid_B::B
end

cross(A::DA, B::DB) where DA<:GridSpace{d1} where DB<:CartesianSpace{d2} where d1 where d2 = ProductSpace{d1+d2, DA, DB}(A,B)

ndims(p::P) where P<:ProductSpace = ndims(p.grid_A) + ndims(p.grid_B)

abstract type Function{d} end

abstract type Linear end

struct MLinear<:Linear end

struct DFun{Dom, Gar, Itp}
    domain::Dom
    values::Gar
    itp::Itp
end


DFun(domain, values) = DFun(domain, values, MLinear())

# Gar<:CartesianGrid
# 1d
function (f::DFun{A,B,I})(i, x::Float64)  where A where B<:GArray{G,V} where V where I<:Linear where G<:PGrid{g1, g2}  where g1<:SSGrid where g2<:CGrid
    ranges = f.values.grid.g2.ranges
    data = f.values[i,:] # this probably allocates
    interp(ranges, data, x)
end

function (f::DFun{A,B,I})(x::Float64)  where A where B<:GArray{G,V} where V where I<:Linear where G<:CGrid 
    ranges = f.values.grid.ranges
    data = f.values.data
    interp(ranges, data, x)
end

ndims(df::DFun) = ndims(df.domain)


# functions defined by hand

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

function discretize(dom::ProductSpace{d,D1,D2}; n=ntuple(u->15,d)) where d where D1<:GridSpace where  D2<:CartesianSpace
    grid_a = SSGrid(dom.grid_A.points)
    g_b = dom.grid_B
    grid_b = CGrid(
        tuple( 
            ( (g_b.min[i], g_b.max[i], n[i]) for i=1:ndims(g_b) )...
        )
    )
    return grid_a × grid_b
end
