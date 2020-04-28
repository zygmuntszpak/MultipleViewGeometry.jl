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




abstract type ProjectiveEntity end
abstract type AbstractContext end
abstract type TwoView <: AbstractContext end

struct HomographyEstimationTask <: AbstractEstimationTask end
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
include("model/plane.jl")
include("model/geometry.jl")
include("model/world.jl")
include("model/projection.jl")
include("model/carrier.jl")
include("model/homography_matrix.jl")
include("model/fundamental_matrix.jl")
include("fit/fit_homography.jl")
include("fit/fit_fundamental_matrix.jl")
include("view/visualize_properties.jl")
include("context/aquire_context.jl")
include("context/visualize_context.jl")



# allotment.jl
export IntervalAllotment

# camera.jl
export  AbstractCamera,
        AbstractCameraModel,
        Pinhole,
        AbstractIntrinsicParameters,
        AbstractExtrinsicParameters,
        IntrinsicParameters,
        ExtrinsicParameters,
        extrinsics,
        intrinsics,
        Camera,
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
export fit_homography

# fit_fundamental_matrix.jl
export fit_fundamental_matrix


end # module
