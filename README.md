AlgebraicNumbers.jl
------

[![Build Status](https://travis-ci.org/anj1/AlgebraicNumbers.jl.svg?branch=master)](https://travis-ci.org/anj1/AlgebraicNumbers.jl)
[![Coverage Status](https://coveralls.io/repos/github/anj1/AlgebraicNumbers.jl/badge.svg?branch=master)](https://coveralls.io/github/anj1/AlgebraicNumbers.jl?branch=master)

#### `sqrt(2)^2 == 2`

This package provides a way of dealing with real and complex numbers exactly and with infinite precision. To see how this works, it's useful to compare with familiar number types like integers and so on. Integer arithmetic (e.g. 2+2==4) is exact but is limited to the operations +, -, and \*. While adding, subtracting, or multiplying two integers always produces another integer, that's not always true with division. With [*rational* numbers](http://docs.julialang.org/en/release-0.4/manual/complex-and-rational-numbers/#rational-numbers), division is included as well. Since many numerical operations can be reduced to sequences of the four elementary operations, this allows a wider range of exact arithmetic to be carried out. *Algebraic* numbers take this further, including not only the four elementary operations, but also *root-taking* operations, for example sqrt() and cbrt(). More generally, the *n*th root of an algebraic number `x` can be taken with:

```julia
root(x, n)
```

And this will be represented exactly. For instance, you can see for yourself that:

```julia
# sqrt(x) is defined as root(x, 2)
sqrt(AlgebraicNumber(2))^2 == 2
```

And this is true for any integer:
```julia
# x = any integer
sqrt(AlgebraicNumber(x))^2 == x
```

Here, `AlgebraicNumber` is just a constructor that takes a number (either an integer or a rational number) and produces an algebraic number.

Note that if you display an algebraic number, you might get something like:
```julia
julia> AlgebraicNumber(1)
≈1.0 + 0.0im
```

That is, something that looks like an approximate complex number, not an exact number. This is *only* the library's way of *displaying* algebraic numbers, and it's simply because in general it is impossible to represent an algebraic number exactly in decimal notation no matter how many digits you display! Internally, algebraic numbers are represented exactly, but they are not represented using decimal or floating-point representation (more on internal representation below).

Indeed, you can do arithmetic on algebraic numbers and all results will be represented exactly:

```julia
sqrt2 = sqrt(AlgebraicNumber(2))
sqrt3 = sqrt(AlgebraicNumber(3))
sqrt6 = sqrt2*sqrt3
# a simple example
assert(sqrt6 == sqrt(AlgebraicNumber(6)))

# slightly more complicated
x = 1 + sqrt6
assert((x - 1)^2 == 6)

# even more complicated
assert(sqrt6 == sqrt(x^2 - 2*sqrt6 - 1))

# and here's another one
y = sqrt(x)
assert((y^2 - 1)^2 == 6)
```

Even *more* generally, arbitrary root-taking operations are possible. That is, you can represent the root of any polynomial (with integer, rational, or algebraic coefficients) as an algebraic number, even if that root doesn't have a representation in terms of a sequence of +, -, /, *, and root-taking operations.

#### Internal implementation

Computer algebra systems (CASes) also allow you to represent algebraic numbers, but the method they use is somewhat different. In CAS systems, numbers are typically represented using the expressions used to generate them. So `sqrt(2)` would be literally represented as `sqrt(2)`. Thus `^2` and `sqrt` cancel out to give `2`. That approach is flexible but it has a fairly large computational cost. The way algebraic numbers are represented here is different - they are represented as discrete roots of minimal polynomials. This approach is a bit more limiting (for example, the `exp()` of an algebraic number is not necessarily an algebraic number) but it is more computationally efficient and allows doing things like equality testing very rapidly and in a way that is always guaranteed to give the correct result, no matter how complicated the algebraic number is. This is something that CAS systems often cannot do.

The tradeoff in using the minimal polynomial representation is that operations like addition and multiplication become non-trivial to compute, since we need to compute a new minimal polynomial, and this involves computation of [resultants](http://specfun.inria.fr/bostan/publications/BoFlSaSc06.pdf) and polynomial factoring. The code for computing resultants has been written in pure julia (in `newton.jl`) and the polynomial factorization is done using the FLINT library, wrapped with the excellent [Nemo.jl](https://github.com/wbhart/Nemo.jl) package. If you are just using this package, though, you usually do not need to worry about any of this.

See [this blog post](https://pseudoprofound.wordpress.com/2016/07/09/some-fun-with-algebraic-numbers/) for some more description and some neat examples.

#### Extra functions

There are a few extra utility functions. For example, `exp_alg(x)` returns exp(iπx), which, assuming x is a rational number, is algebraic. For example:

```julia
# calculate exp(im*pi*2/3) as an algebraic number
x = exp_alg(2//3)
assert(x == sqrt(AlgebraicNumber(-3))/2 - AlgebraicNumber(1)/2)
```

Similarly, `cos_alg(x)` and `sin_alg(x)` return the cosine and sine of πx, which is algebraic if x is rational. These numbers are known as the [trigonometric](https://en.wikipedia.org/wiki/Trigonometric_number) [numbers](https://en.wikipedia.org/wiki/Trigonometric_constants_expressed_in_real_radicals#2.25.C2.B0:_regular_octacontagon_.2880-sided_polygon.29):

```julia
# An example trigonometric number
x = sin_alg(1//8)
y = sqrt(2 - sqrt(AlgebraicNumber(2)))/2
assert(x == y)

# Another example
x = cos_alg(2//5)
y = (sqrt(AlgebraicNumber(5))-1)/4
assert(x == y)
```
