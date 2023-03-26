using NoLib

import Dolo

nolibmodel = include("models/rbc.jl")


dolomodel = Dolo.Model("/home/pablo/.julia/dev/Dolo/examples/models/rbc_mc.yaml")

@time NoLib.time_iteration(nolibmodel; improve=false, verbose=true)

@time Dolo.time_iteration(dolomodel)

module Temp

    import NoLib
    import NoLib: transition, arbitrage
    import Dolo
    using StaticArrays
    using NoLib: CartesianSpace, GridSpace, ×, SSGrid, CGrid, LVectorLike
    using LabelledArrays

    struct DoModel{A,B,C,DM,N} <: NoLib.AModel{A,B,C,N}
        calibration::A
        domain::B
        transition::C
        source::DM
    end


    struct DoDModel{A,B,C,D,DM,N} <: NoLib.ADModel{A,B,C,D,N}
        calibration::A
        domain::B
        grid::C
        transition::D
        source::DM
    end

    

    function DoModel(filename)

        source = Dolo.yaml_import(filename)


        m = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:exogenous])...) 
        s = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:states])...) 
        x = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:controls])...) 
        p = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:parameters])...) 

        calibration = (;m, s, x, p,)

        p1,q1=size(source.exogenous.values)
        p2,q2=size(source.exogenous.transitions)

        Q = SMatrix{p1,q1}(source.exogenous.values)
        points = [Q[i,:] for i=1:size(Q,1)]

        P = SMatrix{p2,q2}(source.exogenous.transitions)

        (;states, min, max) = source.domain

        aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))

        cspace = CartesianSpace(; aargs... )

        domain = GridSpace( points  ) × cspace

        name = try
            Symbol(source.data[:name].value)
        catch e
            :anonymous
        end

        return DoModel{ typeof(calibration), typeof(domain), typeof(P), typeof(source), Val(name)}(calibration, domain, P, source)
        # transition
    end

    function transition(model::DoDModel, m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
        
        v1 = SVector(m...)
        v2  =SVector(s...) 
        v3 =  SVector(x...)
        v4 = SVector(M...)
        v5 = SVector(p...)

        S = Dolo.transition(model.source, SVector(m...), SVector(s...), SVector(x...), SVector(M...), SVector(p...))

        return LVectorLike(model.calibration.s, S)

    end
    
    function arbitrage(model::DoDModel, m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)

        r = Dolo.arbitrage(model.source, SVector(m...), SVector(s...), SVector(x...), SVector(M...), SVector(S...), SVector(X...), SVector(p...))

        return LVectorLike(model.calibration.x, r)
        
    end


    function DoDModel(filename)

        source = Dolo.yaml_import(filename)


        m = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:exogenous])...) 
        s = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:states])...) 
        x = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:controls])...) 
        p = SLVector(; (v=>source.calibration.flat[v] for v in source.symbols[:parameters])...) 

        calibration = (;m, s, x, p,)

        p1,q1=size(source.exogenous.values)
        p2,q2=size(source.exogenous.transitions)

        Q = SMatrix{p1,q1}(source.exogenous.values)
        points = [Q[i,:] for i=1:size(Q,1)]

        P = SMatrix{p2,q2}(source.exogenous.transitions)

        (;states, min, max) = source.domain

        aargs = (states[i]=>(min[i],max[i]) for i in eachindex(states))

        gpsce = GridSpace( points  )
        cspace = CartesianSpace(; aargs... )

        domain = GridSpace( points  ) × cspace

        name = try
            Symbol(source.data[:name].value)
        catch e
            :anonymous
        end


        # dp, pg = Dolo.discretize(dolomodel)
        egrid = Dolo.discretize(source.domain)

        exo = SSGrid( [Q[i,:] for i=1:size(Q,1)] )
        
        args = ( tuple( ((a,b,c)  for (a,b,c) in zip(egrid.min, egrid.max, egrid.n))... ) ) 

        endo = CGrid( args )

        grid = exo × endo


        return DoDModel{ typeof(calibration), typeof(domain), typeof(grid), typeof(P), typeof(source), Val(name)}(calibration, domain, grid, P, source)
        # transition
    end

end

import Main.Temp
import Main.Temp: transition, arbitrage

@time domodel = Main.Temp.DoDModel("/home/pablo/.julia/dev/Dolo/examples/models/rbc_mc.yaml");

(;m,s,x,p) = domodel.calibration

NoLib.transition(domodel, m,s,x,m,p)

import ForwardDiff

g(u) = NoLib.transition(domodel, m,s,u,m,p)
g(x)
ForwardDiff.jacobian(
    g,
    x
)




f(u)= NoLib.arbitrage(domodel, m,s,u,m,s,x,p)
f(x)
ForwardDiff.jacobian(
    f,
    x
)

NoLib.time_iteration(domodel)