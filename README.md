AlgebraicNumbers.jl
------

#### `sqrt(2)^2 == 2`

This package provides a way of dealing with real numbers exactly and with infinite precision. To see how this works, it's useful to compare with familiar number types like integers and so on. Integer arithmetic (e.g. 2+2==4) is exact but is limited to the operations +, -, and \*. Adding, subtracting, or multiplying two integers produces another integer. With [*rational* numbers](http://docs.julialang.org/en/release-0.4/manual/complex-and-rational-numbers/#rational-numbers), division is included as well. That is, if you add, subtract, muliply, or divide two rational numbers, you get another rational number. This is very useful for exact arithmetic. *Algebraic* numbers take this further, including not only the four elementary operations, but also *root-taking* operations, for example sqrt() and cbrt(). More generally, the *n*th root of an algebraic number `x` can be taken with:

```julia
root(x, n)
```

And this will be represented exactly. For instance, you can see for yourself that:

```julia
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

Computer algebra systems (CASes) also allow you to represent numbers in this way, but the method they use is somewhat different. In CAS systems, numbers are typically represented using the expressions used to generate them. So `sqrt(2)` would be literally represented as `sqrt(2)`. Thus `^2` and `sqrt` cancel out to give `2`. That approach is flexible but it has a fairly large computational cost. The way algebraic numbers are represented here is different - they are represented as discrete roots of minimal polynomials. This approach is a bit more limiting (for example, the `exp()` of an algebraic number is not necessarily an algebraic number) but it is more computationally efficient and allows doing things like equality testing very rapidly and in a way that is always guaranteed to give the correct result, no matter how complicated the algebraic number is. This is something that CAS systems often cannot do.

The tradeoff in using the minimal polynomial representation is that operations like addition and multiplication become non-trivial to compute, since we need to compute a new minimal polynomial, and this involves computation of [resultants](http://specfun.inria.fr/bostan/publications/BoFlSaSc06.pdf) and polynomial factoring. The code for computing resultants has been written in pure julia (in `newton.jl`) and the polynomial factorization is done using the FLINT library, wrapped with the excellent [Nemo.jl](https://github.com/wbhart/Nemo.jl) package. If you are just using this package, though, you usually do not need to worry about any of this.

See [this blog post](https://pseudoprofound.wordpress.com/2016/07/09/some-fun-with-algebraic-numbers/) for some more description and some neat examples.
