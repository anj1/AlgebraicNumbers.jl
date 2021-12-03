struct Multiple{Irrational,Q<:Real} <: AbstractIrrational
    c::Q 
end

function Multiple{R}(n::Q) where R where Q <: Real
    Multiple{R, Q}(n)
end

function printr(m::Multiple{R, Q}) where R where Q
    print(R)
end 

function show(io::IO, m::Multiple{R, Q}) where R where Q <: Integer
    show(io, m.c)
	show(io, R)
end

function show(io::IO, m::Multiple{R, Q}) where R where Q <: Rational
    show(io, numerator(m.c))
    print(io, R)
    print(io, "//")
    show(io, denominator(m.c))
end

function show(io::IO, m::Multiple{R, Q}) where R where Q
    print(io, "(")
    show(io, m.c)
    print(io, ")")
	show(io, R)
end

AbstractFloat(x::Multiple{R,Q}) where {R, Q} = (float(R)*float(x.c))::AbstractFloat
function (::Type{T})(x::Multiple{R,Q}) where T <: Number where R where Q
    P = promote_type(T, Q)
    convert(T, convert(P, R)*convert(P,x.c))::T
end
