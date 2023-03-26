using NoLib

import Dolo

nolibmodel = include("models/rbc.jl")


dolomodel = Dolo.Model("/home/pablo/.julia/dev/Dolo/examples/models/rbc_mc.yaml")

@time NoLib.time_iteration(nolibmodel; improve=false, verbose=true)

@time Dolo.time_iteration(dolomodel)