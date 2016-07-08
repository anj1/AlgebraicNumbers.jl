module AlgebraicNumbers

export AlgebraicNumber
export *,+,-,/,^,==,inv

include("algebraic.jl")
#include("elementary.jl")
#include("compose.jl")
include("newton.jl")

end