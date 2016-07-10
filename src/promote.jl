# promotions from integer and rational types to algebraic number types
import Base.promote

# Algebraic number from integer
AlgebraicNumber{T<:Integer}(x::T) = AlgebraicNumber([-x,one(T)], Complex{BigFloat}(x))
# Algebraic number from rational
AlgebraicNumber{T<:Integer}(x::Rational{T}) = AlgebraicNumber([-num(x), den(x)], Complex{BigFloat}(x))

# promotions
promote{T<:Integer}(x::T, y::AlgebraicNumber) = (AlgebraicNumber(x),y)
promote{T<:Integer}(x::AlgebraicNumber,y::T) = (x,AlgebraicNumber(y))
promote{T<:Integer}(x::Rational{T}, y::AlgebraicNumber) = (AlgebraicNumber(x),y)
promote{T<:Integer}(x::AlgebraicNumber,y::Rational{T}) = (x,AlgebraicNumber(y))