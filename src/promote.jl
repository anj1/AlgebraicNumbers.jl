# conversions and promotions from integer and rational types to algebraic number types
import Base.convert 
import Base.promote_rule

# Algebraic number from integer
#convert{T<:Integer}(::Type{AlgebraicNumber}, x::T) = 
#    AlgebraicNumber(BigInt[-x,one(T)], Complex{BigFloat}(x))
# Algebraic number from rational
#convert{T<:Integer}(::Type{AlgebraicNumber}, x::Rational{T}) =
#    AlgebraicNumber(BigInt[-num(x), den(x)], Complex{BigFloat}(x))

# promotions
promote_rule(x::Type{T},          y::Type{AlgebraicNumber{S,F}}) where {T<:Integer,S,F} = AlgebraicNumber
promote_rule(x::Type{Rational{T}},y::Type{AlgebraicNumber{S,F}}) where {T<:Integer,S,F} = AlgebraicNumber

# conversions back
function convert(::Type{Int64},an::AlgebraicNumber)
	c = an.coeff
	if length(c)==2 && abs(c[2])==1
		return convert(Int64, -c[1]*c[2])
	else
		throw(InexactError())
	end
end