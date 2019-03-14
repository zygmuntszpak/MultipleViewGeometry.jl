__precompile__()
module MultipleViewGeometry

using StaticArrays
using GeometryTypes

# Types exported from `types.jl`
export HomogeneousPoint, ProjectiveEntity, FundamentalMatrix, ProjectionMatrix
export EssentialMatrix, HomographyMatrix, HomographyMatrices
export TotalViews, TwoViews
export ProjectionMatrices, FactorisationAlgorithm, Chojnacki, LatentVariables
export CameraModel, Pinhole, CanonicalLens
export EstimationAlgorithm, DirectLinearTransform, Taubin, FundamentalNumericalScheme
export BundleAdjustment
export CostFunction, ApproximateMaximumLikelihood, AML
export HomogeneousCoordinates
export CoordinateSystemTransformation, CanonicalToHartley, HartleyToCanonical
export CovarianceMatrices
export Point2D, Point2DH, Point3D, Point3DH
export HessianApproximation, CanonicalApproximation, CovarianceEstimationScheme
export NoiseModel, GaussianNoise
export RasterSystem, CartesianSystem, OpticalSystem

# Aliases exported from math_aliases.jl
export ⊗, ∑, √

# Functions exported from `operators.jl`.
export 𝑛, ∂𝑛, smallest_eigenpair,vec2antisym, hom⁻¹, hom, ∂hom⁻¹

# Functions exported from `hartley_transformation.jl`.
export hartley_normalization, hartley_normalization!, hartley_transformation

# Functions exported from `transform.jl`.
export transform

# Functions exported from `camera.jl`.
export ascertain_pose, obtain_intrinsics

# Functions exported from `moments_fundamentalmatrix.jl`
export moments

# Functions exported from `estimate_twoview.jl`
export estimate

# Functions exported from `construct_essentialmatrix.jl`
export construct

# Functions exported from `construct_fundamentalmatrix.jl`
export construct

# Functions exported from `construct_projectionmatrix.jl`
export construct

# Functions exported from `project.jl`
export project, epipole

# Functions exported from `move.jl`
export rotx, roty, rotz, rotxyz, rodrigues2matrix
export rotate!, translate!, relocate!

# Functions exported from `cost_functions.jl`
export cost, X, covariance_matrix, covariance_matrix_debug

# Functions exported from "synthetic_data.jl"
export generate_planar_points

# Functions exported from `draw.jl`
#export draw!, EpipolarLineGraphic, LineSegment3D, PlaneSegment3D, Camera3D
#export WorldCoordinateSystem3D
export draw!


# Functions exported from `constraints.jl`
export satisfy, EpipolarConstraint, Constraint

# Functions exported from `triangulation.jl`
export triangulate

# Functions exported from `noise.jl`
export perturb

include("math_aliases/ModuleMathAliases.jl")
include("types/ModuleTypes.jl")
include("operators/ModuleOperators.jl")
#include("rotation/ModuleRotation.jl")
include("synthetic_data/ModuleSyntheticData.jl")
include("camera/ModuleCamera.jl")
include("move/ModuleMove.jl")
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
#using .ModuleRotation
using .ModuleSyntheticData
using .ModuleCamera
using .ModuleMove
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
