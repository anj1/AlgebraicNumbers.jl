module AlgebraicNumbers

using Nemo
import PolynomialRoots

export AlgebraicNumber
export *,+,-,/,^,root,==,inv
export sqrt,cbrt
export exp_alg,cos_alg,sin_alg
export pow2


include("algebraic.jl")
include("promote.jl")
include("newton.jl")

end
