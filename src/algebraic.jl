# Exact representation of algebraic numbers
# (Numbers that are roots of polynomials with integer coefficients)
# And also arithmetic on algebraic numbers,
# including +, -, *, /, and radicals.

using Nemo
import PolynomialRoots.roots

# see: http://nemocas.org/nemo-0.3.pdf

# An algebraic number,
# consisting of the minimal polynomial of the number,
# an arbitrary-precision approximation of the number,
# and prec which specifies the minimal distance between roots of p
# TODO: apprx has to be complex.
type AlgebraicNumber
	coeff::Vector{Int64}
	apprx::Complex{BigFloat}
	prec::BigFloat
end
# algebraic number from just poly and approximation.
# computes precision and simplifies as well.
function AlgebraicNumber(coeff::Vector{Int64}, apprx::Complex{BigFloat})
	an = AlgebraicNumber(coeff, apprx, BigFloat(0.0))
	calc_precision!(an)
	return simplify(an)
end
# Algebraic number from integer
AlgebraicNumber(x::Integer) = AlgebraicNumber([-x,1], Complex{BigFloat}(x))
# Algebraic number from rational
AlgebraicNumber{T<:Integer}(x::Rational{T}) = AlgebraicNumber([-num(x), den(x)], Complex{BigFloat}(x))


function poly_from_coeff(a)
	R,x=PolynomialRing(ZZ,"x")
	sum([a[i]*x^(i-1) for i=1:length(a)])
end

import Base.show
# TODO: only show up to precision
function show(io::IO, an::AlgebraicNumber)
	print(io,"≈")
	#ndigits = max(10, round(Int,ceil(convert(Float64,log(an.prec)/log(10)))))
	show(io,convert(Complex{Float64},an.apprx))
	#print(io,"...")
end

#get_coeffs(p::Nemo.fmpz_poly) = pointer_to_array(convert(Ptr{Int64}, p.coeffs), (p.length,))
get_coeffs(p::Nemo.fmpz_poly) = Int64[Nemo.coeff(p,i) for i=0:degree(p)]
prec_roots(a::Vector{Int64}) = PolynomialRoots.roots(convert(Array{BigFloat},a))
# TODO: make sure roots returns distinct roots

# compute smallest distance between all pairs of elements in x
function min_pairwise_dist(x)
	biginf = convert(BigFloat,Inf)
	n = length(x)
	pdists = [i < j ? abs(x[i]-x[j]) : biginf for i=1:n,j=1:n]
	return minimum(pdists)
end

# # find all roots of a polynomial, with enough precision such that
# # all distinct roots have distinct floating-point representations
# function roots_minprecision(polynomial, options...)
# end

# Given an algebraic number, find minimum precision required
# to specify it among roots of an.p
# TODO: handle case of repeated roots precisely
function calc_precision!(an::AlgebraicNumber)
	# first, find all roots of p
	rts = prec_roots(an.coeff)

	# first, trivial case
	if length(rts)==1
		an.prec = convert(BigFloat, Inf)
		return
	end

	# find minimum pairwise distance between roots;
	# multiply by 0.5 safety factor
	dist = 0.5*min_pairwise_dist(rts)
	an.prec = dist
	return
end


# simplify an algebraic number by reducing p to the minimal polynomial.
# This assumes that calc_precision! has already been called.
function simplify(an::AlgebraicNumber)
	# for all factors of an.p, find the one that matches our roots
	fctrs = keys(Nemo.factor(poly_from_coeff(an.coeff)))

	# first, trivial case
	if length(fctrs)==1
		return AlgebraicNumber(get_coeffs(first(fctrs)),an.apprx,an.prec)
	end
	# case where more than one factor exists
	for fctr in fctrs
		coeff = get_coeffs(fctr)
		# TODO: instead of computing roots, substitute apprx and find closest to zero.
		fctr_rts = prec_roots(coeff)
		dists = Bool[abs(an.apprx - r) < an.prec for r in fctr_rts]
		if any(dists)
			# we've found our winner!
			return AlgebraicNumber(coeff,an.apprx,an.prec)
		end
	end
	# should _not_ be at this point!
	throw(ArgumentError("Precision error!"))
end

import Base.==
import Base.inv
import Base.^
import Base.*
import Base.+
==(an1::AlgebraicNumber,an2::AlgebraicNumber) = (an1.coeff==an2.coeff) && abs(an1.apprx-an2.apprx)<min(an1.prec,an2.prec)
inv(an::AlgebraicNumber) = AlgebraicNumber(reverse(an.coeff), inv(an.apprx))

# interleave each elemnet of a with n zeros
interleave(a,n) =  vec(vcat(a',zeros(Int64,n,length(a))))
function nthroot(an::AlgebraicNumber,n::Int64)
	# # if power is zero, return 1
	# # TODO: handle case where an.apprx=0.0
	# if n==0
	# 	R,x=PolynomialRing(ZZ,"x")
	# 	return AlgebraicNumber(x-1,BigFloat(1.0),BigFloat(1.0))
	# end
	# # if power is negative, first invert.
	# if n < 0
	# 	an = inv(an)
	# 	n = -n
	# end
	if n==0
		throw(ArgumentError("n must be nonzero"))
	end
	if n==1
		return an
	end
	# TODO: negative case
	# TODO: quickly calculate precision
	return AlgebraicNumber(interleave(an.coeff, n-1), an.apprx^(1/n))
end

function pow2(an::AlgebraicNumber)
	cfs = an.coeff 
	cfs2 = [iseven(i) ? -cfs[i] : cfs[i] for i=1:length(cfs)]
	pp = poly_from_coeff(cfs)*poly_from_coeff(cfs2)
	p2 = get_coeffs(pp)[1:2:end]
	return AlgebraicNumber(p2, an.apprx*an.apprx)
end


# partially simplify a polynomial b
# eliminating repeated factors
reduce_repeated_factors(p::Nemo.fmpz_poly) = prod(keys(Nemo.factor(p)))

# multiplication
function *(an1::AlgebraicNumber,an2::AlgebraicNumber)
	# Create polynomial ring over x and z
	R, z = PolynomialRing(ZZ, "z")
	S, x = PolynomialRing(R,  "x")

	# make sure an1=/=an2;
	# otherwise just use pow2
	if an1==an2
		return pow2(an1)
	end 

	# p1 is simply the same as an1.p
	p1 = (z^0)*sum([an1.coeff[i]*(x^(i-1)) for i=1:length(an1.coeff)])
	# p2 prime is an2.p[z/x]*x^degree(an2.p)
	d = length(an2.coeff)
	p2 = sum([an2.coeff[i]*(z^(i-1))*(x^(d-i+1)) for i=1:length(an2.coeff)])
	
	#@show p1, p2
	p = reduce_repeated_factors(resultant(p1,p2))
	#g=gcd(p2,p1)
	#@show "hi"
	#@show g
	#@show resultant(p1,p2)
	return AlgebraicNumber(get_coeffs(p),an1.apprx*an2.apprx)
end
function +(an1::AlgebraicNumber,an2::AlgebraicNumber)
	# Create polynomial ring over x and z
	R, z = PolynomialRing(ZZ, "z")
	S, x = PolynomialRing(R,  "x")
	
	# p1 is simply the same as an1.p
	p1 = sum([an1.coeff[i]*(z^0)*(x^(i-1)) for i=1:length(an1.coeff)])
	# p2 prime is an2.p[z-x]
	p2 = sum([an2.coeff[i]*(z-x)^(i-1) for i=1:length(an2.coeff)])

	p = reduce_repeated_factors(resultant(p1,p2))
	return AlgebraicNumber(get_coeffs(p),an1.apprx+an2.apprx)
end

# take roots of a polynomial,
# and return them as algebraic numbers
function alg_roots(coeff::Vector{Integer})
end

confirm_algnumber(b) = sum(b.coeff .* [b.apprx^(i-1) for i=1:length(b.coeff)])

function test1(n)
	coeff = rand(1:10,n+1)
	a = AlgebraicNumber(coeff, BigFloat(0.0), BigFloat(0.0))
	a.apprx = roots(a.coeff)[rand(1:n)]
	calc_precision!(a)
	a = simplify(a)
	@show a.coeff, convert(Complex{Float64},a.apprx), a.prec
	b = nthroot(a,2)
	#b = a*a
	@show b.coeff, convert(Complex{Float64},b.apprx), b.prec
	#c = nthroot(b,2)
	c = b*b
	#c = pow2(b)
	@show c.coeff, convert(Complex{Float64},c.apprx), c.prec
	@show roots(a.coeff) 
	@show roots(b.coeff) 
	@show roots(c.coeff) 	
end

function test2(n)
	coeff = rand(1:10,n+1)
	a = AlgebraicNumber(coeff, BigFloat(0.0), BigFloat(0.0))
	a.apprx = roots(a.coeff)[rand(1:n)]
	calc_precision!(a)

	coeff = rand(1:10,n+1)
	b = AlgebraicNumber(coeff, BigFloat(0.0), BigFloat(0.0))
	b.apprx = roots(b.coeff)[rand(1:n)]
	calc_precision!(b)	

	c = a*b
	abs(roots(c.coeff) .- c.apprx)
end

# sqrt2 = nthroot(AlgebraicNumber(2),2)
# an.p = (x^2-2)*(x^2-3)
# calc_precision!(an)
# an = simplify(an)
# @show an.p  (should be x^2-2)


# sqrt2 = nthroot(AlgebraicNumber(2),2)
# sqrt3 = nthroot(AlgebraicNumber(3),2)
# sqrt6=sqrt2*sqrt3
# sqrt6_ = nthroot(AlgebraicNumber(6),2)
# calc_precision!(sqrt6_)
# @show sqrt6 == sqrt6_

# this little example has got me stumped.
# an=nthroot(nthroot(AlgebraicNumber(3),2) + AlgebraicNumber(-1),2)
# b = an*an
# b.coeff and b.apprx don't match up!