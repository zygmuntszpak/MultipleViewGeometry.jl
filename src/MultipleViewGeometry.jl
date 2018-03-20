__precompile__()

module MultipleViewGeometry

using Compat

include("types.jl")
include("math_aliases.jl")

# Types exported from `types.jl`
export HomogeneousPoint

# Functions exported from `operators.jl`.
export 𝑛

# Functions exported from `hartley_transformation.jl`.
export hartley_normalization, hartley_transformation

include("operators.jl")
include("data_normalization/hartley_transformation.jl")

# package code goes here

end # module
