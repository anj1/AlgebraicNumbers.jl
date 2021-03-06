using Test
using AlgebraicNumbers 

function test1(n)
	coeff = rand(1:10,n+1)
	a = AlgebraicNumber(coeff, Complex{BigFloat}(0.0), BigFloat(0.0))
	a.apprx = prec_roots(a.coeff)[rand(1:n)]
	calc_precision!(a)
	a = simplify(a)
	@show a.coeff, convert(Complex{Float64},a.apprx), a.prec
	b = root(a,2)
	#b = a*a
	@show b.coeff, convert(Complex{Float64},b.apprx), b.prec
	#c = root(b,2)
	c = b*b
	#c = pow2(b)
	@show c.coeff, convert(Complex{Float64},c.apprx), c.prec
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

function test3()
	a = sqrt(AlgebraicNumber(2))
	b = sqrt(AlgebraicNumber(3))
	axb = a*b
	@show axb, axb.coeff
	apb = a+b 
	@show apb, apb.coeff
end 

function test4()
	# simple test
	sqrt2 = sqrt(AlgebraicNumber(2))
	@test sqrt2^2 == 2

	# Golden ratio
	ϕ = 1//2 + sqrt(AlgebraicNumber(5)/4)

	# As we all know, this has the property that:
	@test 1+1/ϕ == ϕ
end

function plastic_constant_test()
	# see http://mathworld.wolfram.com/PlasticConstant.html
	a = sqrt(AlgebraicNumber(69))
	n = cbrt(9-a) + cbrt(9+a)
	p = n*inv(cbrt(AlgebraicNumber(18)))

	@test p-1==1/(p^4)
	@test p+1==p^3
end

function test5()
	# just an answer to a stackoverflow Q.
	# http://math.stackexchange.com/questions/422233/how-to-find-a-minimal-polynomial-field-theory
	n = sqrt(AlgebraicNumber(9*5))-sqrt(AlgebraicNumber(4*7))+sqrt(AlgebraicNumber(35))
	d = 1-sqrt(AlgebraicNumber(5))+sqrt(AlgebraicNumber(7))
	α=n/d
	@test α.coeff == BigInt[3596, 2312, -280, -156, 19]
end

function test_abs()
	ii = sqrt(AlgebraicNumber(-1))
	@test conj(ii) == -ii
	@test abs(ii) == 1
	@test abs(AlgebraicNumber(-7//8))==7//8
end

function test_real_imag()
	a = root(AlgebraicNumber(-1),5)
	alg_im = sqrt(AlgebraicNumber(-1))
	@test real(a) + alg_im*imag(a) == a

	x = rand(0:20)//10
	@test cos_alg(x - 1//2) == sin_alg(x)
	x = 2//1
	@test cos_alg(x - 1//2) == sin_alg(x)
end

# TODO: add some more tests
function test_pow2()
	a = AlgebraicNumber(3//2)
	@test pow2(a) == AlgebraicNumber(9//4)
end

function test_show()
	a = IOBuffer()
	show(a, sqrt(AlgebraicNumber(-1))+1)

	@test convert(UTF8String, takebuf_array(a)) == "≈1.0 + 1.0im"
end

test4()
test5()
plastic_constant_test()
test_abs()
test_real_imag()
test_pow2()
#test_show()

# sqrt2 = root(AlgebraicNumber(2),2)
# an.p = (x^2-2)*(x^2-3)
# calc_precision!(an)
# an = simplify(an)
# @show an.p  (should be x^2-2)


# sqrt2 = root(AlgebraicNumber(2),2)
# sqrt3 = root(AlgebraicNumber(3),2)
# sqrt6=sqrt2*sqrt3
# sqrt6_ = root(AlgebraicNumber(6),2)
# calc_precision!(sqrt6_)
# @show sqrt6 == sqrt6_

# this little example has got me stumped.
# an=root(root(AlgebraicNumber(3),2) + AlgebraicNumber(-1),2)
# b = an*an
# b.coeff and b.apprx don't match up!