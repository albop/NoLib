include("models/neoclassical.jl")

using NoLib: DFun

i0 = 3
s0_ = [NoLib.enum(model.grid)...][i0]
s0 = [NoLib.enum(model.grid)...][i0]
m0 = s0[2:end]


x = GVector(model.grid, [Iterators.repeated(SVector(model.calibration.x), length(model.grid))...])

### TODO: interpolation not defined for pgrid

φ = DFun(model.grid, x)

φ.itp.θ

x0 = x[3]

S = [NoLib.τ(model, s0, x0)...][1][2]

r = NoLib.F(model, s0, x0, φ)



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
