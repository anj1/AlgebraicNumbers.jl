module AlgebraicNumbers

using Nemo
import PolynomialRoots

export AlgebraicNumber
export *,+,-,/,^,root,==,inv
export sqrt,cbrt
export exp_alg,cos_alg,sin_alg
export log_alg,acos_alg,asin_alg
export pow2


include("algebraic.jl")
include("trig.jl")
include("promote.jl")
include("newton.jl")
include("inv_totient.jl")

end
