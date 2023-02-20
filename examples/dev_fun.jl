using NoLib
const NL = NoLib
import NoLib: ×
using StaticArrays


# dom = NL.CartesianSpace(;
#     x=(-10.0, 10.0),
#     y=(0.0, Inf)
# )

dom = NL.CartesianSpace(;x=(-10.0,10.0))
@assert isbits(dom)
grid = NL.discretize(dom)
@assert isbits(grid)


# interpolate scalar values
vals = NL.GVector(grid, [e[1]^2 for e in grid])

φ_l = NL.DFun(dom, vals; interp_mode=:linear)
v = SVector(0.1)
@time φ_l(v)

φ_c = NL.DFun(dom, vals; interp_mode=:cubic)
@time φ_c(v)

using Plots
plot(φ_l)
plot!(φ_c)

vv = [v for i=1:1000]

@time φ_l(vv)
@time φ_c(vv)



# interpolate vector values
vals = NL.GVector(grid, [SVector(e[1]^2) for e in grid])

φ_l = NL.DFun(dom, vals; interp_mode=:linear)
v = SVector(0.1)
@time φ_l(v)

φ_c = NL.DFun(dom, vals; interp_mode=:cubic)
@time φ_c(v)

vv = [v for i=1:1000]

@time φ_l(vv)
@time φ_l.(vv)
@time φ_c(vv)
@time φ_c.(vv)




### 2 dimension

dom = NL.CartesianSpace(;x=(-15.0,15.0), y=(-10.0,10.0))
@assert isbits(dom)
grid = NL.discretize(dom)
@assert isbits(grid)

vals = NL.GVector(grid, [e[1]^2+0.5e[2]^2 for e in grid])


vv = SVector(0.1, 0.2)
vvec = [vv for i=1:1000]

φ_l = NL.DFun(dom, vals; interp_mode=:linear)
@time φ_l(vv)

φ_c = NL.DFun(dom, vals ; interp_mode=:cubic)
@time φ_c(vv)


using Plots
xr = range(grid.ranges[1][1:2]..., 1000)
yr = range(grid.ranges[2][1:2]..., 1000)
z = Surface((x,y)->φ_l(SVector(x,y))-x^2-y^2, xr, yr)
surface(xr,yr,z)


φ_l(SVector(0.1, 0.5))

# @time φ_c.(vvec)
# @time φ_c(vvec)



# mixed grid


dom = NL.GridSpace(
    SVector(
        SVector(-0.1, 0.1),
        SVector(-0.1, -0.1)
    )
) × NL.CartesianSpace(;x=(-10.0,10.0))
@assert isbits(dom)
grid = NL.discretize(dom)
@assert isbits(grid)
vals = NL.GVector(grid, [ SVector(e[1]^2) for e in grid])


φ_l = NL.DFun(dom, vals; interp_mode=:linear)
v = SVector(0.1)
@time φ_l(v)

φ_c = NL.DFun(dom, vals; interp_mode=:cubic)
@time φ_c(v)

vv = [v for i=1:1000]

@time φ_l(vv)
@time φ_c(vv)

















@time φ_c(vv)