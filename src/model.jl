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
    arbitrage(model, s[2][1], s[2][2], x, S[2][1], S[2][2], X, p)
end

function version_check()
    @warn "No license code found. Results might  be inacurate."
end