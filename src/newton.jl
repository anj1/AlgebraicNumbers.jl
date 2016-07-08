using Nemo 

# derivative of polynomial
derivative(c::Vector) = c[2:end] .* (1:length(c)-1)

function polyinv{T}(coeffs::Vector{T}, n)
	R, x = Nemo.PowerSeriesRing(QQ, n, "x")
	a = R(map(QQ, coeffs), length(coeffs), n)
	ai = inv(a)
	return Nemo.fmpq[coeff(ai,i) for i=0:n-1]
end

# compute newton power series of polynomial given with coefficients coeff
function to_newton(coeffs::Vector{BigInt},n)
	# first, make monic.
	coeffs = coeffs//coeffs[end]

	d = length(coeffs)-1
	a_cfs = reverse(derivative(coeffs))
	b_cfs = reverse(coeffs)
	
	# initialize power series polynomials
	R, x = Nemo.PolynomialRing(QQ, "x")
	a = R(map(QQ, a_cfs))
	b = R(map(QQ, b_cfs))
	b0 = R(polyinv(b_cfs, n))

	c  = truncate(a*b0, d)

	r = R()
	x_power = R(1)
	x_d = x^d

	l = round(Int64,floor(n/d))
	for j = 0 : l
		r += c*x_power
		x_power *= x_d
		#c = -truncate(shift_right(b*c,d)*b0,d)
		c = mullow(shift_right(b*c,d),b0,d)
	end
	return Rational{BigInt}[Rational(coeff(r,i)) for i=0:n]
end


# tr: traces i.e. newton series
function from_newton{T}(tr::Vector{T})
	c = T[]
	for k = 1 : length(tr)-1
		push!(c, -dot(tr[2:(k+1)], vcat(reverse(c),1))/k)
	end
	# strip leading zeros
	#c = c[1:findlast(c)]
	return vcat(reverse(c),1)
end 


# composed product of two polynomials, given as coeffs p and q
function composed_product(p::Vector{BigInt},q::Vector{BigInt})
	# compute newton series
	n = (length(p)-1)*(length(q)-1)+1
	a = to_newton(p,n)
	b = to_newton(q,n)

	# multiply newton series and invert
	pq = from_newton(a.*b)

	# convert to integer and return
	return map(num, pq*lcm(map(den, pq)))
end