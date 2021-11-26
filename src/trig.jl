# compute exp(pi*i*q),
# which is algebraic if q is rational.
function exp_alg(q::Rational)
	# first, obtain polynomial
	p = interleave(BigInt[-1,1], 2*denominator(q)-1)
	# now, select root.
	apprx = exp(im*BigFloat(pi)*q)
	# Finally, return minimal polynomial w.r.t. that root
	return AlgebraicNumber(p,apprx)
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