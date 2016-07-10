# Exact representation of algebraic numbers
# (Numbers that are roots of polynomials with integer coefficients)
# And also arithmetic on algebraic numbers,
# including +, -, *, /, and radicals.

using Nemo
import PolynomialRoots
import PolynomialRoots:roots

# see: http://nemocas.org/nemo-0.4.pdf

# An algebraic number,
# consisting of the minimal polynomial of the number,
# an arbitrary-precision approximation of the number,
# and prec which specifies the minimal distance between roots of p
# TODO: apprx has to be complex.
type AlgebraicNumber{T<:Integer,F<:AbstractFloat} <: Number
	coeff::Vector{T}
	apprx::Complex{F}
	prec::F
end

# algebraic number from just poly and approximation.
# computes precision and simplifies as well.
function AlgebraicNumber{T,F}(coeff::Vector{T}, apprx::Complex{F})
	an = AlgebraicNumber{T,F}(coeff, apprx, zero(F))
	calc_precision!(an)
	return simplify(an)
end

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
get_coeffs(p::Nemo.fmpz_poly) = [BigInt(Nemo.coeff(p,i)) for i=0:Nemo.degree(p)]
prec_roots{T<:Integer}(a::Vector{T}) = PolynomialRoots.roots(convert(Array{BigFloat},a))
# TODO: make sure roots returns distinct roots

# Given an algebraic number, find minimum precision required
# to specify it among roots of an.p
# TODO: handle case of repeated roots precisely
function calc_precision!(an::AlgebraicNumber)
	# compute smallest distance between all pairs of elements in x
	function min_pairwise_dist(x)
		biginf = convert(BigFloat,Inf)
		n = length(x)
		pdists = [i < j ? abs(x[i]-x[j]) : biginf for i=1:n,j=1:n]
		return minimum(pdists)
	end

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
	R, x = PolynomialRing(ZZ, "x")
	p = R(map(ZZ, an.coeff))
	fctr_dict = Nemo.factor(p)
	fctrs = keys(fctr_dict)

	# first, trivial case
	if length(fctrs)==1
		# irreducible case
		if first(values(fctr_dict))==1
			return AlgebraicNumber(get_coeffs(first(fctrs)),an.apprx,an.prec)
		end
		# reducible case
		an = AlgebraicNumber(get_coeffs(first(fctrs)),an.apprx,an.prec)
		calc_precision!(an)  # re-calculate precision
		return an
	end
	# case where more than one factor exists
	mindists = [minimum(abs(an.apprx .- prec_roots(get_coeffs(fctr)))) for fctr in fctrs]
	(newprec, i) = findmin(mindists)
	fctr = collect(fctrs)[i]
	an = AlgebraicNumber(get_coeffs(fctr),an.apprx,newprec)
	return an 
end

import Base.==
import Base.inv
import Base.^
import Base.*
import Base.+
import Base.-
import Base./
function ==(an1::AlgebraicNumber,an2::AlgebraicNumber)
	cf1 = an1.coeff
	cf2 = an2.coeff
	(cf1/cf1[end])==(cf2/cf2[end]) || return false 
	calc_precision!(an1)
	calc_precision!(an2)
	return abs(an1.apprx-an2.apprx)<min(an1.prec,an2.prec)
end

inv(an::AlgebraicNumber) = AlgebraicNumber(reverse(an.coeff), inv(an.apprx))

# interleave each elemnet of a with n zeros
interleave(a,n) =  vec(vcat(a',zeros(Int64,n,length(a))))
function root(an::AlgebraicNumber,n::Int64)
	if n==0
		throw(ArgumentError("n must be nonzero"))
	end
	if n==1
		return an
	end
	if n < 0
		an = inv(an)
		n = -n
	end
	# TODO: quickly calculate precision
	return AlgebraicNumber(interleave(an.coeff, n-1), an.apprx^(1/n))
end

import Base.sqrt 
import Base.cbrt
sqrt(an::AlgebraicNumber) = root(an,2)
cbrt(an::AlgebraicNumber) = root(an,3)

# TODO: special, more efficient cases for ^2 and ^3
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
	p = composed_product(an1.coeff, an2.coeff)
	return AlgebraicNumber(p, an1.apprx * an2.apprx)
end

function +(an1::AlgebraicNumber,an2::AlgebraicNumber)
	p = composed_sum(an1.coeff, an2.coeff)
	return AlgebraicNumber(p, an1.apprx + an2.apprx)
end

function -(an1::AlgebraicNumber)
	cfs = an1.coeff
	for i=1:2:length(cfs)
		cfs[i]=-cfs[i]
	end
	return AlgebraicNumber(cfs, -an1.apprx, an1.prec)
end

-(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1+(-an2)
/(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1*(inv(an2))

# take roots of a polynomial,
# and return them as algebraic numbers
function alg_roots(coeff::Vector{Integer})
	#TODO
end

confirm_algnumber(b) = sum(b.coeff .* [b.apprx^(i-1) for i=1:length(b.coeff)])