# compute exp(pi*i*q),
# which is algebraic if q is rational.
function exp_alg(q::Rational)
	if numerator(q)==0
		return AlgebraicNumber(1)
	end

	# first, obtain minimal polynomial
	s, x = PolynomialRing(AlgebraicNumbers.Nemo.ZZ, "x")
	poly = cyclotomic(2*denominator(q), x)
	coeffs = [convert(BigInt, coeff(poly, i)) for i=0: 2*denominator(q)]

	# now, select root.
	apprx = exp(im*BigFloat(pi)*q)

	# Finally, return polynomial w.r.t. that root
	return AlgebraicNumber(coeffs, apprx)
end

cos_alg(q::Rational) = real(exp_alg(q))
sin_alg(q::Rational) = imag(exp_alg(q))

# checks if a polynomial is equal to any of the 
# cyclotomic polynomials n in a list of candidates
function is_cyclotomic(poly, candidates, x)
	for n in candidates
		if poly == cyclotomic(n, x)
			return (true, n)
		end 
	end
	return (false, 0)
end 

# compute log(a)/(pi*i),
# which is rational if a is a root of unity.
# If a is not a root of unity, returns Nothing
function log_alg(a::AlgebraicNumber)
	s, x = PolynomialRing(Nemo.ZZ, "x")

	deg = length(a.coeff)-1
	poly = s(map(Nemo.ZZ, a.coeff))

	(is_cycl, denom) = is_cyclotomic(poly, inv_totient(deg), x)

	if is_cycl
		# TODO: should this be BigInt?
		num = round(Int, denom*imag(log(a.apprx)/pi))
		return num//denom
	else 
		return Nothing 
	end
end

function acos_alg(x::AlgebraicNumber)
	# TODO: check if the number can be made a root of unity.

	# First, make the number a root of unity.
	y = sqrt(1 - x^2)
	z = x + sqrt(AlgebraicNumber(-1))*y

	# Now take log 
	log_alg(z)
end

function asin_alg(x::AlgebraicNumber)
	# TODO: check if the number can be made a root of unity.

	# First, make the number a root of unity.
	y = sqrt(1 - x^2)
	z = y + sqrt(AlgebraicNumber(-1))*x

	# Now take log 
	log_alg(z)
end