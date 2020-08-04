# _Fast_ composed sums and composed products of polynomials,
# using the algorithm described in:
# "Fast computation of special resultants"
# by Bostan, Flajolet, Salvy, and Schost

# derivative of polynomial
derivative(c::Vector) = c[2:end] .* (1:length(c)-1)

function polyinv(coeffs::Vector, n)
	R, x = Nemo.PowerSeriesRing(Nemo.FlintQQ, n, "x")
	a = R(map(Nemo.FlintQQ, coeffs), length(coeffs), n, 0)
	ai = inv(a)
	return Nemo.fmpq[coeff(ai,i) for i=0:n-1]
end

# compute newton power series of polynomial given with coefficients coeff,
# in base field R,x.
# See fig.1 in reference
function to_newton(coeffs::Vector{BigInt},n,R,x)
	# first, make monic.
	coeffs = coeffs//coeffs[end]

	d = length(coeffs)-1
	a_cfs = reverse(derivative(coeffs))
	b_cfs = reverse(coeffs)

	# initialize power series polynomials
	a = R(map(Nemo.FlintQQ, a_cfs))
	b = R(map(Nemo.FlintQQ, b_cfs))
	b0 = R(polyinv(b_cfs, n))

	c  = truncate(a*b0, d)

	r = R()
	x_power = R(1)
	x_d = x^d

	l = round(Int64,floor(n/d))
	for j = 0 : l
		r += c*x_power
		x_power *= x_d
		c = -mullow(shift_right(b*c,d),b0,d)
	end
	return r
end

to_array(p) = Rational{BigInt}[Rational(coeff(p,i)) for i=0:Nemo.degree(p)]

# tr: traces i.e. newton series
# This algorithm is based on the Leverrier-Faddeev algorithm
# see: http://math.stackexchange.com/questions/405822/what-is-the-fastest-way-to-find-the-characteristic-polynomial-of-a-matrix
function from_newton(tr::Vector{T}) where {T<:Number}
	# special case
	if tr==[1]
		return T[0,1]
	end
	n = length(tr)
	c = Array{T}(UndefInitializer(),n)
	c[end] = one(T)
	for k = 1 : n-1
		next_c = -sum(tr[2:(k+1)].*c[end-k+1:end])/k
		c[end-k] = next_c
	end
	return c
end

# Hadamard (element-wise) product of two polynomials
function hadm(p,q,R)
	n = min(Nemo.degree(p),Nemo.degree(q))
	R([Nemo.coeff(p,i)*Nemo.coeff(q,i) for i=0:n])
end

# composed product of two polynomials, given as coeffs p and q
function composed_product(p::Vector{BigInt},q::Vector{BigInt})
	# compute newton series
	n = (length(p)-1)*(length(q)-1)+1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	a = to_newton(p,n,R,x)
	b = to_newton(q,n,R,x)

	# multiply newton series and invert
	pq = from_newton(to_array(hadm(a,b,R)))

	# convert to integer and return
	return map(numerator, pq*lcm(map(denominator, pq)))
end

# composed sum of two polynomials, given as coeffs p and q
function composed_sum(p::Vector{BigInt},q::Vector{BigInt})
	# compute newton series
	n = (length(p)-1)*(length(q)-1)+1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	a = to_newton(p,n,R,x)
	b = to_newton(q,n,R,x)

	# exp series
	ee  = R([Nemo.FlintQQ(1//factorial(BigInt(i))) for i=0:n])
	eei = R([Nemo.FlintQQ(   factorial(BigInt(i))) for i=0:n])

	# multiply newton series and invert
	m = truncate(hadm(a,ee,R)*hadm(b,ee,R),n+1)
	pq = from_newton(to_array(hadm(m,eei,R)))

	# convert to integer and return
	return map(numerator, pq*lcm(map(denominator, pq)))
end
