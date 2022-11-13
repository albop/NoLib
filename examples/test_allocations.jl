include("neoclassical_model.jl")


φ = GVector(model.grid, [Iterators.repeated(SVector(model.calibration.x), length(model.grid))...])
r = deepcopy(φ)


function check_indexing(model, φ)
    # sum( model.grid[i][2] for i=1:length(model.grid) )
    # n = length(model.grid)
    n = length(model.grid)
    n
end

@time check_indexing(model, φ);

function check_cover()
    a = SVector(1.0, 2.0, 3.0, 4.0)
    m = SVector(0.0, 0.0)
    c = NoLib.cover(m,a)
    sum(c)
end

@time check_cover()

function compute_transition(model, φ)
    i = 1
    s = ((1,1),(model.grid[1]))
    r = sum(s for (w,(i,s)) in NoLib.τ(model, s, x))
    return sum(r)
end


function compute_res!(r, model, φ)
    NoLib.F!(r, model, φ, φ)
end

function compute_res_2!(r, model, φ)
    for (i,s) in NoLib.enum(model.grid)
        x = φ[i...]
        r[i...] = NoLib.F(model, (i,s), x, φ)
    end
end


@time compute_transition(model, φ);




compute_res!(r, model, φ)
compute_res_2!(r, model, φ)

@assert (@allocated compute_res!(r, model, φ)) == 0
@assert (@allocated compute_res_2!(r, model, φ)) == 0
