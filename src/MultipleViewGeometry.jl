__precompile__()

module MultipleViewGeometry

using Compat

include("types.jl")
include("math_aliases.jl")

# Types exported from `types.jl`
export HomogeneousPoint, ProjectiveEntity, FundamentalMatrix

# Functions exported from `operators.jl`.
export 𝑛, smallest_eigenpair

# Functions exported from `hartley_transformation.jl`.
export hartley_normalization, hartley_transformation

# Functions exported from `moments_fundamentalmatrix.jl`
export moments

# Functions exported from `estimate_twoview.jl`
export estimate

include("operators.jl")
include("data_normalization/hartley_transformation.jl")
include("twoview_estimation/moments_fundamentalmatrix.jl")
include("twoview_estimation/estimate_twoview.jl")

# package code goes here

end # module
