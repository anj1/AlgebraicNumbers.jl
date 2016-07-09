module AlgebraicNumbers

export AlgebraicNumber
export *,+,-,/,^,root,==,inv
export sqrt,cbrt

include("algebraic.jl")
include("newton.jl")

end