__precompile__()

module MultipleViewGeometry

using Compat

include("types.jl")
include("math_aliases.jl")

# Types exported from `types.jl`
export HomogeneousPoint, ProjectiveEntity, FundamentalMatrix, ProjectionMatrix
export CameraModel, Pinhole, CanonicalLens

# Functions exported from `operators.jl`.
export 𝑛, smallest_eigenpair,vec2antisym

# Functions exported from `hartley_transformation.jl`.
export hartley_normalization, hartley_transformation

# Functions exported from `moments_fundamentalmatrix.jl`
export moments

# Functions exported from `estimate_twoview.jl`
export estimate

# Functions exported from `construct_fundamentalmatrix.jl`
export construct

# Functions exported from `construct_projectionmatrix.jl`
export construct

# Functions exported from `project.jl`
export project

# Functions exported from `rotations.jl`
export rotx, roty, rotz, rotxyz, rodrigues2matrix

include("operators.jl")
include("rotation/rotations.jl")
include("data_normalization/hartley_transformation.jl")
include("camera/construct_projectionmatrix.jl")
include("projection/project.jl")
include("twoview_estimation/moments_fundamentalmatrix.jl")
include("twoview_estimation/estimate_twoview.jl")
include("twoview_estimation/construct_fundamentalmatrix.jl")

# package code goes here

end # module
