# Exact representation of algebraic numbers
# (Numbers that are roots of polynomials with integer coefficients)
# And also arithmetic on algebraic numbers,
# including +, -, *, /, and radicals.

import PolynomialRoots

import Base.zero,Base.one
import Base.+,Base.-,Base.*,Base./,Base.inv
import Base.abs,Base.conj
import Base.real,Base.imag
import Base.==

# see: http://nemocas.org/nemo-0.4.pdf

# An algebraic number,
# consisting of the minimal polynomial of the number,
# an arbitrary-precision approximation of the number,
# and prec which specifies the minimal distance between roots of p
# TODO: apprx has to be complex.
struct AlgebraicNumber{T<:Integer,F<:AbstractFloat} <: Number
	coeff::Vector{T}
	apprx::Complex{F}
	prec::F
end

# algebraic number from just poly and approximation.
# computes precision and simplifies as well.
function AlgebraicNumber(coeff::Vector{T}, apprx::Complex{F}) where {T<:Integer,F<:AbstractFloat}
	an = AlgebraicNumber{T,F}(coeff, apprx, calc_precision(coeff, apprx))
	return simplify(an)
end

# AlgebraicNumber from any integer type
AlgebraicNumber(x::T) where {T<:Integer} =
    AlgebraicNumber(BigInt[-x,one(T)], Complex{BigFloat}(x))

# AlgebraicNumber from rationails
AlgebraicNumber(x::Rational) =
    AlgebraicNumber(BigInt[-numerator(x), denominator(x)], Complex{BigFloat}(x))

AlgebraicNumber(x::Complex) =
    AlgebraicNumber(real(x)) + AlgebraicNumber(imag(x))*root(AlgebraicNumber(-1),2)


function poly_from_coeff(a)
	R,x=polynomial_ring(Nemo.FlintZZ,"x")
	sum([a[i]*x^(i-1) for i=1:length(a)])
end

function is_displayed_exactly(an)
	io = IOBuffer()
	show(io,convert(Complex{Float64},an.apprx))
	displ = String(take!(io))
	from_displ = AlgebraicNumber(Complex{Rational{BigInt}}(parse(Complex{BigFloat}, displ)))
	from_displ==an, displ
end

import Base.show
# TODO: only show up to precision
function show(io::IO, an::AlgebraicNumber)
	display_exact, displ = is_displayed_exactly(an)
	display_exact || print(io, "≈")
	print(io, displ)
end

#get_coeffs(p::Nemo.ZZPolyRingElem) = pointer_to_array(convert(Ptr{Int64}, p.coeffs), (p.length,))
get_coeffs(p::Nemo.ZZPolyRingElem) = [BigInt(Nemo.coeff(p,i)) for i=0:Nemo.degree(p)]
prec_roots(a::Vector{T}) where {T<:Integer} = PolynomialRoots.roots(convert(Array{BigFloat},a))
# TODO: make sure roots returns distinct roots

# Given an algebraic number, find minimum precision required
# to specify it among roots of an.p
# TODO: handle case of repeated roots precisely
function calc_precision(coeff::Vector{T}, apprx::Complex{F}) where {T<:Integer,F<:AbstractFloat}
	# compute smallest distance between all pairs of elements in x
	function min_pairwise_dist(x)
		biginf = convert(F,Inf)
		n = length(x)
		if n<=1
			return Inf
		else
			pdists = [i < j ? abs(x[i]-x[j]) : biginf for i=1:n,j=1:n]
			return minimum(pdists)
		end
	end

	# first, find all roots of p
	rts = prec_roots(coeff)

	# first, trivial case
	if length(rts)==1
		return convert(F, Inf)
	end

	# find minimum pairwise distance between roots;
	# multiply by 0.5 safety factor
	return 0.5*min_pairwise_dist(convert(Vector{Complex{F}},rts))
end


# simplify an algebraic number by reducing p to the minimal polynomial.
# This assumes that calc_precision! has already been called.
function simplify(an::AlgebraicNumber)
	# for all factors of an.p, find the one that matches our roots

	# If linear polynomial, then already irreducible.
	if length(an.coeff)<=2
		return an
	end

	# Otherwise, factor out.
	R, x = polynomial_ring(Nemo.FlintZZ, "x")
	p = R(map(Nemo.FlintZZ, an.coeff))
	fctr_dict = Nemo.factor(p)
	#fctrs = keys(fctr_dict)
	fctrs = [p for (p,e) in fctr_dict]

	# first, trivial case
	if length(fctrs)==1
		# irreducible case
		if first(values(fctr_dict))==1
			return AlgebraicNumber(get_coeffs(first(fctrs)),an.apprx,an.prec)
		end
		# reducible case
		coeffs1 = get_coeffs(first(fctrs))
		apprx1  = an.apprx
		return AlgebraicNumber(coeffs1,apprx1,calc_precision(coeffs1,apprx1))
	end
	# case where more than one factor exists
	mindists = [minimum(abs.(an.apprx .- prec_roots(get_coeffs(fctr)))) for fctr in fctrs]
	(newprec, i) = findmin(mindists)
	fctr = collect(fctrs)[i]
	an = AlgebraicNumber(get_coeffs(fctr),an.apprx,newprec)
	return an
end

function ==(an1::AlgebraicNumber,an2::AlgebraicNumber)
	cf1 = an1.coeff
	cf2 = an2.coeff
	(cf1/cf1[end])==(cf2/cf2[end]) || return false
	prec1 = calc_precision(an1.coeff, an1.apprx)
	prec2 = calc_precision(an2.coeff, an2.apprx)
	return abs(an1.apprx-an2.apprx)<min(prec1,prec2)
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
	# first check if it is already in the form of a square root.
	if all(cfs[2:2:end] .== 0)
		pp_cfs = cfs
	else
		cfs2 = [iseven(i) ? -cfs[i] : cfs[i] for i=1:length(cfs)]
		pp = poly_from_coeff(cfs)*poly_from_coeff(cfs2)
		pp_cfs = get_coeffs(pp)
	end
	p2 = pp_cfs[1:2:end]
	return AlgebraicNumber(p2, an.apprx*an.apprx)
end


# partially simplify a polynomial b
# eliminating repeated factors
reduce_repeated_factors(p::Nemo.ZZPolyRingElem) = prod(keys(Nemo.factor(p)))

# multiplication
function *(an1::AlgebraicNumber,an2::AlgebraicNumber)
	if an1==0 || an2==0
		# TODO: don't handle this explicitly
		return zero(AlgebraicNumber)
	end
	# check if p==q, if then use a more optimized and correct routine
	#if an1.coeff == an2.coeff
	#	return
	#end
	p = composed_product(an1.coeff, an2.coeff)
	return AlgebraicNumber(p, an1.apprx * an2.apprx)
end

function +(an1::AlgebraicNumber,an2::AlgebraicNumber)
	p = composed_sum(an1.coeff, an2.coeff)
	return AlgebraicNumber(p, an1.apprx + an2.apprx)
end

function -(an1::AlgebraicNumber)
	cfs = copy(an1.coeff)
	for i=1:2:length(cfs)
		cfs[i]=-cfs[i]
	end
	return AlgebraicNumber(cfs, -an1.apprx, an1.prec)
end

-(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1+(-an2)
/(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1*(inv(an2))

# the complex conjugate of an algebraic number has the same minimal polynomial
conj(an::AlgebraicNumber) = AlgebraicNumber(an.coeff,conj(an.apprx),an.prec)
abs(an::AlgebraicNumber) = sqrt(an*conj(an))

zero(::Type{AlgebraicNumber}) = AlgebraicNumber(BigInt[0, 1],Complex{BigFloat}(0.0),BigFloat(1.0))
one(::Type{AlgebraicNumber})  = AlgebraicNumber(BigInt[-1,1],Complex{BigFloat}(1.0),BigFloat(1.0))

real(an::AlgebraicNumber) = (an+conj(an))*AlgebraicNumber(BigInt[1,-2], BigFloat(0.5)+0im,BigFloat(0.5))
imag(an::AlgebraicNumber) = (an-conj(an))*AlgebraicNumber(BigInt[1,0,4],BigFloat(-0.5)*im,BigFloat(0.5))

# take roots of a polynomial,
# and return them as algebraic numbers
function alg_roots(coeff::Vector{Integer})
	#TODO
end

confirm_algnumber(b) = sum(b.coeff .* [b.apprx^(i-1) for i=1:length(b.coeff)])

