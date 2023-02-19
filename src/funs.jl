abstract type Domain{d} end
    
struct CartesianDomain{d}
    # n::NTuple{d, Int64}
    min::NTuple{d, Float64}
    max::NTuple{d, Float64}
end
ndims(cd::CartesianDomain{d}) where d = d

struct GridDomain{N,d}
    points::SVector{N,SVector{d,Float64}}
end
ndims(gd::GridDomain{d}) where d = d

GridDomain(v::Vector{SVector{d, Float64}}) where d = GridDomain(SVector(v...))

struct ProductDomain{d,A,B}
    grid_A::A
    grid_B::B
end

cross(A::DA, B::DB) where DA<:GridDomain{d1} where DB<:CartesianDomain{d2} where d1 where d2 = ProductDomain{d1+d2, DA, DB}(A,B)

ndims(p::P) where P<:ProductDomain = ndims(p.grid_A) + ndims(p.grid_B)

abstract type Function{d} end

struct DFun{Dom, Gar}
    domain::Dom
    values::Gar
end


ndims(df::DFun) = ndims(df.domain)

# function ddFun(dom, gar)
#     d = ndims(dom)
#     dt = typeof(dom)
#     gt = typeof(gar)
#     DFun{d, dt, gt}(dom, gar)
# end

function discretize(dom::ProductDomain{d,D1,D2}; n=ntuple(u->15,d)) where d where D1<:GridDomain where  D2<:CartesianDomain
    grid_a = SSGrid(dom.grid_A.points)
    g_b = dom.grid_B
    grid_b = CGrid(
        tuple( 
            ( (g_b.min[i], g_b.max[i], n[i]) for i=1:ndims(g_b) )...
        )
    )
    return grid_a × grid_b
end
