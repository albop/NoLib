using NoLib

include("models/rbc.jl")

@time mem = NoLib.time_iteration_workspace(model);
@time begin NoLib.time_iteration_4(model, mem; verbose=false,T=5) end
@time begin NoLib.time_iteration_4(model, mem; verbose=false, improve=true) end


#mystery: why is time_iteration_1 twice as fast?
# Response: it calls F 170 times instead of 71
@time begin NoLib.time_iteration_1(model; verbose=false) end
@time begin NoLib.time_iteration_4(model; verbose=false) end

#mystery: why is time_iteration_1 twice as fast?
@benchmark begin NoLib.time_iteration_1(model; verbose=false) end
@benchmark begin NoLib.time_iteration_4(model; verbose=false) end


@profview begin NoLib.time_iteration_1(model; verbose=false) end
@profview begin NoLib.time_iteration_4(model; verbose=false) end


@code_warntype NoLib.time_iteration_1(model; verbose=false)
@code_warntype NoLib.time_iteration_4(model, mem; verbose=false)


@time begin NoLib.time_iteration_4(model; verbose=true) end
@time begin NoLib.time_iteration_4(model; verbose=true, improve=true) end




@time mem = NoLib.time_iteration_workspace(model);

@time begin NoLib.time_iteration_4(model, mem; verbose=false) end
