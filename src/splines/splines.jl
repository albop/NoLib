module splines

export filter_coeffs, interpolant_cspline, filter_coeffs
export eval_UC_spline, eval_UC_spline!
export prefilter!

include("csplines.jl")
include("splines_filter.jl")
include("cubic_prefilter.jl")


function interpolant_cspline(a, b, orders, V)

    coefs = filter_coeffs(a, b, orders, V)

    function fun(s::Array{Float64,2})
        return eval_UC_spline(a, b, orders, coefs, s)
    end

    function fun(p::Float64...)
        return fun([p...]')
    end

    return fun

end

abstract type Linear end
abstract type Cubic end

struct MLinear<:Linear end
struct MCubic<:Cubic end

struct CubicInterpolator{G,C} <: Cubic
    grid::G
    θ::C
end


function CubicInterpolator(grid; values=nothing)

    n = [e[3] for e in grid.ranges]
    θ = zeros(eltype(values), (i+2 for i in n)...)
    ci = CubicInterpolator{typeof(grid), typeof(θ)}(grid, θ)
    if !isnothing(values)
        ind = tuple( (2:(e[3]+1) for e in grid.ranges )...)
        ci.θ[ind...] .= values
        splines.prefilter!(ci.θ)
    end
    return ci

end


end
