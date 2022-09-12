include("rbc_def.jl")


sol = NoLib.time_iteration_3(model;verbose=true, improve=false);

sol = NoLib.time_iteration_1(model;verbose=true);


x0 = sol.solution


a = model.grid.g2.ranges[1][1]
b = model.grid.g2.ranges[1][2]
xval = [range( a,b; length=100)...]

yval = [x0(2, SVector(e))[1] for e in xval]



using Plots

(;δ) = model.p



plot(xval, xval*δ)

a,b=([e[2][1] for e in NoLib.iti(model.grid.g2)], [x0[2,j][1] for j=1:length(model.grid.g2)])
scatter!(a,b)

scatter!([model.s[1]],[model.x[1]])
plot!(xval, yval)










using StaticArrays

grid = model.grid

l = [SVector(model.x) for e in 1:length(grid)]

x0 = GArray(
    grid,
    l
)


using NoLib: enum
using StaticArrays
#
m = SVector(grid[1][1:1]...)
s = SVector(grid[1][2:2]...)

a = x0[1]

ss  = ( (1,2), (m,s))

gen = NoLib.τ(model, ss::Tuple, a::SVector)

[gen...]


gen = NoLib.τ_fit(model, ss::Tuple, a::SVector)

[gen...]


μ0 = GArray(grid, ones(length(grid))/length(grid))

@time μ1 = NoLib.G(model, μ0, x0)


@time NoLib.transition_matrix(model, x0)


[NoLib.τ_fit(model,(ss[1][1], ss[2]), a)...]


x0(ss)

NoLib.F(model, ss, x0[1], x0)

NoLib.F(model, x0, x0);

@time A = NoLib.dF_1(model, x0, x0);

@time NoLib.dF_2(model, x0, x0, x0*0.000001);

@time sol = NoLib.time_iteration_1(model;);
@time sol = NoLib.time_iteration_2(model;);



sol = NoLib.time_iteration_3(model;verbose=true, improve=false);

# @time sol = NoLib.time_iteration_3(model;verbose=true, improve=true);


@time NoLib.time_iteration_3(model;verbose=true, improve=true);








using Plots
# @time convert(Matrix, L);

# spy(convert(Matrix, L))

A = convert(Matrix, NoLib.dF_1(model, x0, x0));

B = NoLib.dF_2(model, x0, x0);

⟂(a,b) = min(a,b)
⟂ᶠ(a,b) = min(a,b)

34 - 33 ⟂ 4

⊥