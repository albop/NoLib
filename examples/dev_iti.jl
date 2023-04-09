using NoLib

include("models/rbc.jl")

@time mem = NoLib.time_iteration_workspace(model);
@time begin NoLib.time_iteration(model, mem; verbose=false) end


@time mem = NoLib.time_iteration_workspace(model);
@time begin NoLib.time_iteration(model, mem; verbose=false, T=5) end
@time begin NoLib.time_iteration(model, mem; verbose=true, improve=true) end


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
@time begin NoLib.time_iteration(model, mem; verbose=false, T=10, interp_mode=:cubic) end;
@time begin NoLib.time_iteration(model, mem; verbose=false, improve=true, interp_mode=:cubic, T=10) end;


@time mem = NoLib.time_iteration_workspace(model);

@time begin NoLib.time_iteration_4(model, mem; verbose=false) end
