__precompile__()
module MultipleViewGeometry

using Compat
using StaticArrays



# Types exported from `types.jl`
export HomogeneousPoint, ProjectiveEntity, FundamentalMatrix, ProjectionMatrix
export CameraModel, Pinhole, CanonicalLens
export EstimationAlgorithm, DirectLinearTransform, Taubin, FundamentalNumericalScheme
export CostFunction, ApproximateMaximumLikelihood, AML
export HomogeneousCoordinates
export CoordinateSystemTransformation, CanonicalToHartley, HartleyToCanonical
export CovarianceMatrices
export Point2DH, Point3DH
export HessianApproximation, CanonicalApproximation, CovarianceEstimationScheme
export NoiseModel, GaussianNoise

# Aliases exported from math_aliases.jl
export ⊗, ∑, √

# Functions exported from `operators.jl`.
export 𝑛, ∂𝑛, smallest_eigenpair,vec2antisym

# Functions exported from `hartley_transformation.jl`.
export hartley_normalization, hartley_normalization!, hartley_transformation

# Functions exported from `transform.jl`.
export transform

# Functions exported from `moments_fundamentalmatrix.jl`
export moments

# Functions exported from `estimate_twoview.jl`
export estimate

# Functions exported from `construct_fundamentalmatrix.jl`
export construct

# Functions exported from `construct_projectionmatrix.jl`
export construct

# Functions exported from `project.jl`
export project, epipole

# Functions exported from `rotations.jl`
export rotx, roty, rotz, rotxyz, rodrigues2matrix

# Functions exported from `cost_functions.jl`
export cost, X, covariance_matrix, covariance_matrix_debug

# Functions exported from `draw.jl`
export draw!, EpipolarLineGraphic, LineSegment3D, PlaneSegment3D, Camera3D
export WorldCoordinateSystem3D

# Functions exported from `constraints.jl`
export satisfy, EpipolarConstraint, Constraint

# Functions exported from `triangulation.jl`
export triangulate

# Functions exported from `noise.jl`
export perturb

include("math_aliases/ModuleMathAliases.jl")
include("types/ModuleTypes.jl")
include("operators/ModuleOperators.jl")
include("rotation/ModuleRotation.jl")
include("data_normalization/ModuleDataNormalization.jl")
include("transform/ModuleTransform.jl")
include("projection/ModuleProjection.jl")
include("carriers/ModuleCarriers.jl")
include("moments/ModuleMoments.jl")
include("cost_function/ModuleCostFunction.jl")
include("estimate/ModuleEstimation.jl")
include("construct/ModuleConstruct.jl")
include("draw/ModuleDraw.jl")
include("constraints/ModuleConstraints.jl")
include("triangulation/ModuleTriangulation.jl")
include("noise/ModuleNoise.jl")

using .ModuleMathAliases
using .ModuleTypes
using .ModuleOperators
using .ModuleRotation
using .ModuleDataNormalization
using .ModuleTransform
using .ModuleProjection
using .ModuleCarriers
using .ModuleEstimation
using .ModuleMoments
using .ModuleCostFunction
using .ModuleConstruct
using .ModuleDraw
using .ModuleConstraints
using .ModuleTriangulation
using .ModuleNoise


# package code goes here

end # module
