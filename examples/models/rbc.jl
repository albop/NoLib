using NoLib

using StaticArrays
using LabelledArrays
using NoLib: SGrid, CGrid, PGrid, GArray, GVector, enum, SSGrid
import NoLib: transition, arbitrage
import NoLib: ×, DModel

## Define Model

model = let
	
	β = 0.9
	σ = 5
	η = 1
	δ = 0.025
	α = 0.33
	ρ = 0.8
	zbar = 0
	σ_z = 0.016
	n = 0.33
	z = zbar
	rk = 1/β - 1+δ
	k = n/(rk/α)^(1/(1-α))
	w = (1-α)*exp(z)*(k/n)^α
	y = exp(z)*k^α*n^(1-α)
	i = δ*k
	c = y - i
	χ =  w/c^σ/n^η
	
	p = (; β, σ, η, χ, δ, α, ρ, zbar, σ_z)
	
	m = SLVector( (; z))
    s = SLVector( (; k))
    x = SLVector( (; n, i))
	
	P = @SMatrix [0.4 0.6; 0.6 0.4]
	Q = @SMatrix [-0.01; 0.01]
	
	exo = SSGrid( [Q[i,:] for i=1:size(Q,1)] )
    endo = CGrid( ((0.5*k, 1.5*k, 100),) )
    grid = exo × endo
	
	DModel(
		(;m, s, x, p,),
		grid,
		P
	)
	
end

function intermediate(model::typeof(model), m::SLArray, s::SLArray, x::SLArray, p)
	y = exp(m.z)*(s.k^p.α)*(x.n^(1-p.α))
	w = (1-p.α)*y/x.n
	rk = p.α*y/s.k
	c = y - x.i
	return SLVector( (; y, c, rk, w))
end

function transition(model::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    K = s.k*(1-p.δ) + x.i
    return SLVector( (;K) )
end
	
function arbitrage(model::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
	y = intermediate(model, m, s, x, p)
	Y = intermediate(model, M, S, X, p)
	res_1 = p.χ*(x.n^p.η)*(y.c^p.σ) - y.w
	res_2 = (p.β*(y.c/Y.c)^p.σ)*(1 - p.δ + Y.rk) - 1
    return SLVector( (;res_1, res_2) )
end