# Composed products and composed sums of polynomials
# see "Fast computation of special resultants" by Bostan et. al. 2006

using Polynomials
using TaylorSeries

typealias Coeff Float64 #Rational{Int64}

import Base.round 
round{T<:Integer}(::Type{Poly{T}},p::Poly) = Poly(round(T,coeffs(p)))

# reverse poly
polyrev(p::Poly) = Poly(reverse(coeffs(p)))

function polytrunc(p::Poly,l::Int,u::Int)
	cfs = coeffs(p)
	assert(l >= 0)
	low_bnd = l+1
	upp_bnd = min(u,length(cfs))
	return Poly(cfs[low_bnd:upp_bnd])
end
polyceil( p::Poly,u::Int) = polytrunc(p,0,u)
polyfloor(p::Poly,l::Int) = polytrunc(p,l,degree(p)+1)

# hadamard product
import Base:.*
.*(f::Poly,g::Poly) = Poly(coeffs(f).*coeffs(g))

# invert polynomial
# e.g. 1/p
# to precision n
function polyinv(p::Poly,n::Integer)
	#cfs = convert(Array{Rational{Int64}},coeffs(p))
	#tp = Taylor1{Rational{Int64}}(cfs,n) 
	cfs = convert(Array{Coeff},coeffs(p))
	tp = Taylor1{Coeff}(cfs,n) 
	inv_tp = 1/tp
	return Poly(inv_tp.coeffs[1:n])
end
# TODO: see Sieveking–Kung algorithm (Sieveking, 1972; Kung, 1974),

# see fig. 1, pp. 6
function newton(h::Poly,n::Integer)
	d = degree(h)
	a = polyrev(polyder(h))
	b = polyrev(h)
	b0 = polyinv(b,d)
	c = polytrunc(a*b0,0,d)
	l = round(Int64,floor(n/d))
	r = Poly([0])
	x_d = Poly([0,1])^d
	x_power = Poly([1])
	for j = 0 : l
		r += c*x_power
		x_power *= x_d
		c=-polyceil(polyfloor(b*c,d)*b0,d)
	end
	# TODO: change to polyceil
	return polyceil(r,n)
end

# this function is used as a 'gold standard'
# for verifying the above function.
function newton_slow(h::Poly,n::Integer)
	p = polyrev(polyder(h))*polyinv(polyrev(h),n)
	return polyceil(p,n)
end


function poly_divr(p::Poly,n::Integer)
	cfs = convert(Array{Coeff},coeffs(p))
	tp = Taylor1{Coeff}(cfs,n) 
	cfs = convert(Array{Coeff},coeffs(polyder(p)))
	tpp = Taylor1{Coeff}(cfs,n) 
	t = tpp/tp
	return Poly(t.coeffs[1:n])
end

function from_newton(nh::Poly,d)
	x = Poly([0,1])
	#s = div(nh - coeffs(nh)[1], x)
	s = polyfloor(nh,1)
	r = 1 - coeffs(s)[1]*x
	n = 2
	while (n<=d)
		#m_prime = -polyceil(polyder(r)*polyinv(r,2n)+s, 2n-1)
		m_prime = -polyceil(poly_divr(r,2n)+s, 2n-1)
		dmp = degree(m_prime)
		m = 1 + Poly(coeffs(m_prime).*[i==0 ? 0//1 : 1//i for i=0:dmp])
		r = polyceil(r*m, 2n)
		n = 2n 
		@show r
	end
	r = polyceil(r,d+1)
	return polyrev(r)
end

function from_newton_slow(nh::Poly,n)
	s = polyint(-polyfloor(nh,1))
	cfs = convert(Array{Coeff},coeffs(s))
	tp = Taylor1{Coeff}(cfs,n)
	etp = exp(tp)
	Poly(reverse(etp.coeffs[1:n]))
end

# Lemma 3, pp. 10
# TODO: do something about monic conversion
# TODO: use faster from_newton
# TODO: use rational coefficients when possible.
function composed_product(f,g)
	nm = degree(f)*degree(g)
	from_newton_slow(newton(f,nm+1).*newton(g,nm+1),nm+1)
end


### tests
function test_polytrunc()
	p = Poly([0,0,3])
	assert(polyfloor(p,0) == p)
	assert(polyfloor(p,1) == Poly([0,3]))
	assert(polyfloor(p,2) == Poly([3]))
	assert(polyfloor(p,3) == zero(Poly{Int64}))
	assert(polyfloor(p,4) == zero(Poly{Int64}))

	assert(polyceil(p,0) == zero(Poly{Int64}))
	assert(polyceil(p,1) == zero(Poly{Int64}))
	assert(polyceil(p,2) == zero(Poly{Int64}))
	assert(polyceil(p,3) == p)
	assert(polyceil(p,4) == p)
end

# TODO: test from_newton as well
# (once it works, that is!)
# TODO: also handle non-monics!
function test_newton()
	for i = 1 : 100
		degr = rand(1:10)
		cfs = vcat(rand(-10:10,degr),1)
		p = Poly(cfs)
		p2 = round(Poly{Int},from_newton_slow(newton(p,degr+1),degr+1))
		assert(p == p2)
	end
end

# test composed product by producing random polys,
# 