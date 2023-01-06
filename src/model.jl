# functions to help with model definition

struct Model{A,B,C}
    name::Symbol
    calibration::A
    domain::B
    transition::C
end

struct DModel{A,B,C}
    calibration::A
    grid::B
    transition::C
end

function recalibrate()
end

# TODO
import Base.show
Base.show(io::IO, dmodel::DModel) = print(io, "DModel(#$(hash(typeof(dmodel))))")


exo_transition(model::DModel{A,B,C}) where A where B where C = model.transition

import Base: merge
function merge(a::SLArray, b::SLArray)
    # TODO: type is incorrect (should correspond to SVector)
    m = (merge(convert(NamedTuple,a), convert(NamedTuple,b)))
    return SLVector(m)
end

function LVectorLike(m0::SLArray{Tuple{d}, T, 1, d, M}, m) where d where T where M
    tt = eltype(m)
    TT = SLArray{Tuple{d}, tt, 1, d, M}
    return TT(m...)
end

label_GArray(m, g::GArray) = GArray(g.grid, [LVectorLike(m, e) for e in g.data])

function transition(model, m, s, x, M, p)
    m = LVectorLike(model.calibration.m,m)
    s = LVectorLike(model.calibration.s,s)
    x = LVectorLike(model.calibration.x,x)
    M = LVectorLike(model.calibration.m,M)
    S = transition(model,m,s,x,M,p)
    return SVector(S...)
end



function arbitrage(model, m, s, x, M, S, X, p)
    m = LVectorLike(model.calibration.m, m)  # this does not keep the original type
    s = LVectorLike(model.calibration.s, s)
    x = LVectorLike(model.calibration.x, x)
    M = LVectorLike(model.calibration.m, M)
    S = LVectorLike(model.calibration.s, S)
    X = LVectorLike(model.calibration.x, X)
    r = arbitrage(model,m,s,x,M,S,X,p)
    return SVector(r...)
end

function arbitrage(model, s, x, S, X)
    p = model.calibration.p
    
    arbitrage(model, 
        NoLib.split_states(model, s)...,
        x,
        NoLib.split_states(model, S)...,
        X,
        p
    )
end



function split_states(model, s_)

    n_m = length(model.calibration.m)  # this does not keep the original type
    n_s = length(model.calibration.s)

    m = SVector((s_[i] for i=1:n_m)...)
    s = SVector((s_[i] for i=n_m+1:(n_m+n_s))...)

    return (;m,s)

end

## default implementation
function initial_guess(model, m::SLArray, s::SLArray, p)
    model.calibration.x
end

function initial_guess(model, m::SVector, s::SVector, p)

    m = LVectorLike(model.calibration.m, m)  # this does not keep the original type
    s = LVectorLike(model.calibration.s, s)
    
    x = initial_guess(model, m, s, p)
    
    return SVector(x...)
    
end

function initial_guess(model, s::SVector)
    
    p = model.calibration.p

    m_, s_ = split_states(model, s)

    x = initial_guess(model, m_, s_, p)

    return SVector(x...)

end


function initial_guess(model)
    GVector(
        model.grid,
        [initial_guess(model, s) for s in model.grid]
    )
end
