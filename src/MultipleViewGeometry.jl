module MultipleViewGeometry

using Colors
#using GeometryTypes
using GeometryBasics
using InplaceOps
using LinearAlgebra
import Makie
import Makie: AbstractPlotting, AbstractPlotting.Scene, AbstractPlotting.RGBAf0, AbstractPlotting.RGBf0
import Makie: mesh!, scatter!, linesegments!
using OffsetArrays
using Parameters
using PGFPlotsX
using ProjectiveNumericalSchemes
using StaticArrays
using Statistics
using Setfield
using Rotations
using FiniteDiff
using ForwardDiff # TODO add this to test set instead for debugging purposes
using Random


import ProjectiveNumericalSchemes: evaluate_gradient, evaluate_hessian
import ProjectiveNumericalSchemes: HartleyNormalizeDataContext
import ProjectiveNumericalSchemes: construct_carrier, construct_carrier_derivatives
import ProjectiveNumericalSchemes: UncertainObservations
import ProjectiveNumericalSchemes: construct_𝐗
import ProjectiveNumericalSchemes: vec2antisym
import ProjectiveNumericalSchemes: matrices
import ProjectiveNumericalSchemes: hom, hom⁻¹
import ProjectiveNumericalSchemes: ApproximateMaximumLikelihood
import ProjectiveNumericalSchemes: AbstractUncertainObservations
import ProjectiveNumericalSchemes: VectorValuedObjective
import ProjectiveNumericalSchemes: SumOfSquares

const ReprojectionError = SumOfSquares

abstract type ProjectiveEntity end
abstract type AbstractContext end
abstract type TwoView <: AbstractContext end

struct HomographyEstimationTask <: AbstractEstimationTask end
struct CameraCalibrationTask <: AbstractEstimationTask end
struct FundamentalMatrixEstimationTask <: AbstractEstimationTask end
abstract type AbstractFittingMethod end

# const Point = GeometryTypes.Point
# const Vec3 = GeometryTypes.Vec3
# const Point3 = GeometryTypes.Point3
# const Vec = GeometryTypes.Vec

#struct TwoViewExperiment <: TwoView end

include("util.jl")
include("rotations.jl")
include("model/coordinate_system.jl")
include("model/allotment.jl")
include("model/image.jl")
include("model/camera.jl")
include("model/pose.jl")
include("model/plane.jl")
include("model/geometry.jl")
include("model/world.jl")
include("model/projection.jl")
include("model/carrier.jl")
include("model/homography_matrix.jl")
include("model/fundamental_matrix.jl")
include("model/essential_matrix.jl")
include("model/noise.jl")
include("fit/fit_sole_camera_rig.jl")
include("fit/fit_homography.jl")
include("fit/fit_fundamental_matrix.jl")
include("view/visualize_properties.jl")
include("context/aquire_context.jl")
include("context/calibration_context.jl")
include("context/world_system_transformation_context.jl")
include("context/visualize_context.jl")
include("context/triangulate_context.jl")
include("context/rectification_context.jl")


# allotment.jl
export IntervalAllotment

# camera.jl
export  AbstractCamera,
        AbstractCameraModel,
        AbstractLensCameraModel,
        AbstractDistortionModel,
        Pinhole,
        RadialDistortionModel,
        AbstractIntrinsicParameters,
        AbstractExtrinsicParameters,
        IntrinsicParameters,
        ExtrinsicParameters,
        extrinsics,
        intrinsics,
        Camera,
        coefficients,
        distortion,
        model,
        image_type,
        translate,
        relocate

# coordinate_system.jl
export  AbstractCoordinateSystem,
        AbstractPlanarCoordinateSystem,
        PlanarCartesianSystem,
        CartesianSystem,
        RasterSystem,
        OpticalSystem,
        coordinate_system,
        origin,
        basis_vectors

# geometry.jl
export  EuclideanPlane3D,
        AbstractPlane,
        Line3D,
        normal,
        origin,
        distance,
        on_plane,
        on_line,
        on_line

# image.jl
export  AbstractImage,
        AbstractAnalogueImage,
        AbstractDigitalImage,
        AnalogueImage,
        coordinate_system

# plane.jl
export  AbstractPlane,
        Plane,
        PlaneSegment

# projection.jl
export  Projection,
        project,
        back_project,
        ascertain_pose,
        matrix

# homography_matrix.jl
export HomographyMatrix,
       HomographyMatrices,
       matrices

# fundamental_matrix.jl
export FundamentalMatrix

# essential_matrix.jl
export EssentialMatrix

# world.jl
export PrimitiveWorld,
       PlanarWorld

# rotations.jl
export rotxyz,
       rotx,
       roty,
       rotz,
       rodrigues2matrix


# util.jl
export  hom⁻¹,
        hom,
        vec2antisym,
        smallest_eigenpair

# aquire_context.jl
export AquireImage,
       aquire

# visualize_properties.jl
export AbstractVisualProperties,
       MakieVisualProperties,
       PGFPlotsVisualProperties

# visualize_context.jl
export VisualizeWorld

#
export DirectLinearTransform,
       Observations,
       UncertainObservations

# fit_homography.jl
export fit_homography,
       LevenbergMarquardt # TODO move


# fit_fundamental_matrix.jl
export fit_fundamental_matrix

# calibration_context.jl
export CalibrateCamera

# noise.jl
export apply_noise

# pose.jl
export RelativePose,
       CoordinateTransformation,
       rotation,
       translation

# world_system_transformation_context.jl
export WorldSystemTransformation

# triangulate_context.jl
export Triangulate,
       DirectLinearTriangulation

export Rectify,
       FusielloCalibratedRectification

end # module
