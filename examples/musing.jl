include("neoclassical_model.jl")


i0 = 3
s0_ = [NoLib.enum(model.grid)...][i0]
s0 = [NoLib.enum(model.grid)...][i0]
m0 = s0[2:end]


φ = GVector(model.grid, [Iterators.repeated(SVector(model.x), length(model.grid))...])

x0 = φ[3]
S = [NoLib.τ(model, s0, x0)...][1][2]


r = NoLib.dF_1(model, φ, φ)

function check_alloc(model, r, φ)
    NoLib.dF_1!(r, model, φ, φ)
    nothing
end

@time check_alloc(model, r, φ)

r = deepcopy(φ)

function check_alloc(model, φ, r)
    NoLib.dF_2!(r, model, φ, φ, φ)
    nothing
end

@time check_alloc(model,  φ, r)
