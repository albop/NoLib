using NoLib

model = include("models/rbc.jl")

@time mem = NoLib.time_iteration_workspace(model);
@time begin NoLib.time_iteration(model, mem; verbose=false, T=50) end;

@time mem = NoLib.time_iteration_workspace(model);
@time begin NoLib.time_iteration(model, mem; verbose=false, interp_mode=:cubic) end;

@time mem = NoLib.time_iteration_workspace(model);
@time NoLib.time_iteration(model, mem; verbose=true, improve=true) ;


@time begin NoLib.time_iteration(model; verbose=true, improve=true) end


@time begin NoLib.time_iteration(model; verbose=false) end
@time begin NoLib.time_iteration(model; verbose=true, interp_mode=:cubic) end
@time begin NoLib.time_iteration(model; verbose=true, interp_mode=:linear) end



@time begin NoLib.time_iteration(model; verbose=false, interp_mode=:linear) end;
@time begin NoLib.time_iteration(model; verbose=false, interp_mode=:cubic) end;


@time begin NoLib.time_iteration(model; improve=true, verbose=false, interp_mode=:linear) end

@time begin NoLib.time_iteration(model; improve=true, verbose=false, interp_mode=:cubic) end


@time mem = NoLib.time_iteration_workspace(model);
@time begin NoLib.time_iteration(model, mem; verbose=false, T=5) end;
@time begin NoLib.time_iteration(model, mem; verbose=false, improve=true) end;


@time mem = NoLib.time_iteration_workspace(model);
@time sol = begin NoLib.time_iteration(model, mem; verbose=false, T=10, interp_mode=:cubic) end;


@time mem = NoLib.time_iteration_workspace(model);

@time begin NoLib.time_iteration_4(model, mem; verbose=false) end





mem = NoLib.time_iteration_workspace(model, interp_mode=:cubic);
@time sol = begin NoLib.time_iteration(model, mem; verbose=false) end;

mem = NoLib.time_iteration_workspace(model, interp_mode=:linear);
@time sol = begin NoLib.time_iteration(model, mem; verbose=false) end;
    

# @time sol = begin NoLib.time_iteration(model, mem; improve=true, verbose=false) end;

@time mem = NoLib.time_iteration_workspace(model, interp_mode=:cubic);
@time sol = NoLib.time_iteration(model, mem);
    
@time mem = NoLib.time_iteration_workspace(model, interp_mode=:cubic);
@time sol = NoLib.time_iteration(model, mem; improve=true);

@time mem = NoLib.time_iteration_workspace(model, interp_mode=:linear);
@time sol = NoLib.time_iteration(model, mem);
    

@time mem = NoLib.time_iteration_workspace(model, interp_mode=:linear);
@time sol = NoLib.time_iteration(model, mem; improve=true);



x0 = sol.solution
φ = NoLib.DFun(model, x0; interp_mode=:cubic)

typeof(φ)
typeof(φ.itp)

f(φ, x0) = begin NoLib.fit!(φ, x0); nothing end

@time f(φ, x0)

r0 = deepcopy(x0)

@time NoLib.F!(r0, model, x0, φ)

J0 = NoLib.dF_1(model, x0, φ)

@time NoLib.dF_1!(J0, model, x0, φ);

