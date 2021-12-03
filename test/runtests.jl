using Test
using AlgebraicNumbers
import AlgebraicNumbers.inv_totient

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
	# @show axb, axb.coeff
	@test abs(AlgebraicNumbers.confirm_algnumber(axb)) < 1e-10
	apb = a+b
	# @show apb, apb.coeff
	@test abs(AlgebraicNumbers.confirm_algnumber(apb)) < 1e-10
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
	@test String(take!(a)) == "1.0 + 1.0im"

	a = IOBuffer()
	show(a, sqrt(AlgebraicNumber(2)))
	@test String(take!(a))[1] == '≈'
end

function test_log_alg()
	@test log_alg(exp_alg(1//9)) == 1//9

	@test log_alg(AlgebraicNumber(2)) == Nothing
end

function test_trig_alg()
	@test acos_alg(cos_alg(3//7)) == 3//7
	@test asin_alg(sin_alg(3//7)) == 3//7

	@test acos_alg(AlgebraicNumber(1)) == 0//1
	@test asin_alg(AlgebraicNumber(1)) == 1//2

	@test asin_alg(AlgebraicNumber(3//2)) == Nothing
	@test acos_alg(AlgebraicNumber(3//2)) == Nothing
end 

function totient(x::T) where T <: Integer 
	prod([(fac.first^(fac.second-1))*(fac.first-1) for fac in Nemo.factor(x)])
end

# Test the correctness of inv_totient for all totients up to m
function check_inv_totient(m::T) where T <: Integer
    # Lower bound on phi(n)==m,
    # And thus worst-case maximum range we need to consider.
    n = 2*m^2

    tots = [totient(i) for i=1:n]

    for i = 1:m
        _gold = findall(==(i), tots)
        _test = sort(collect(inv_totient(i)))

        length(_gold) == length(_test) || return false

        all(_gold .== _test) || return false
    end 

	return true
end

function test_inv_totient(m::T) where T <: Integer 
	@test check_inv_totient(m)
end 

test3()
test4()
test5()
plastic_constant_test()
test_abs()
test_real_imag()
test_pow2()
test_show()
test_log_alg()
test_trig_alg()

# testcase of issue #5
@test AlgebraicNumber(1)+sqrt(AlgebraicNumber(-1)) != AlgebraicNumber(2)

# sqrt2 = root(AlgebraicNumber(2),2)
# an.p = (x^2-2)*(x^2-3)
# calc_precision!(an)
# an = simplify(an)
# @show an.p  (should be x^2-2)

# test multiplication of square roots
sqrt2 = root(AlgebraicNumber(2),2)
sqrt3 = root(AlgebraicNumber(3),2)
sqrt6=sqrt2*sqrt3
sqrt6_ = root(AlgebraicNumber(6),2)
@test sqrt6 == sqrt6_

#
an=root(root(AlgebraicNumber(3),2) + AlgebraicNumber(-1),2)
b = an*an
@test abs(AlgebraicNumbers.confirm_algnumber(b)) < 1e-10
