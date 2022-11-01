include("neoclassical_model.jl")


φ = GVector(model.grid, [Iterators.repeated(SVector(model.calibration.x), length(model.grid))...])
r = deepcopy(φ)

function compute_res!(r, model, φ)
    NoLib.F!(r, model, φ, φ)
end

function compute_res_2!(r, model, φ)
    for (i,s) in NoLib.enum(model.grid)
        x = φ[i...]
        r[i...] = NoLib.F(model, (i,s), x, φ)
    end
end


compute_res!(r, model, φ)
compute_res_2!(r, model, φ)

@assert (@allocated compute_res!(r, model, φ)) == 0
@assert (@allocated compute_res_2!(r, model, φ)) == 0
