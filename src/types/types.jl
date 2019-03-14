# mutable struct Point2DH <: FieldVector{3,Float64}
#     x::Float64
#     y::Float64
#     h::Float64
# end
#
# mutable struct Point3DH <: FieldVector{4,Float64}
#     x::Float64
#     y::Float64
#     z::Float64
#     h::Float64
# end

#const Point2DH = MMatrix{3,1,Float64,3}
#const Point3DH = MMatrix{4,1,Float64,4}

const Point2DH = SVector{3,Float64}
const Point3DH = SVector{4,Float64}

const Point2D = SVector{2,Float64}
const Point3D = SVector{3,Float64}

# mutable struct Point2DH <: MMatrix{3,1,Float64,3}
#     x::Float64
#     y::Float64
#     h::Float64
# end
#
# mutable struct Point3DH <: MMatrix{4,1,Float64,4}
#     x::Float64
#     y::Float64
#     z::Float64
#     h::Float64
# end


struct HomogeneousPoint{T <: AbstractFloat,N}
    coords::NTuple{N, T}
end

abstract type ProjectiveEntity end

abstract type CameraModel end

abstract type EstimationAlgorithm end

abstract type FactorisationAlgorithm end

abstract type CostFunction end

abstract type CoordinateSystemTransformation end

abstract type CovarianceEstimationScheme end

abstract type NoiseModel end

abstract type TotalViews end

abstract type PlanarCoordinateSystem end


mutable struct FundamentalMatrix <: ProjectiveEntity
end

mutable struct HomographyMatrix <: ProjectiveEntity
end

struct HomographyMatrices <: ProjectiveEntity
end

mutable struct EssentialMatrix <: ProjectiveEntity
end

mutable struct  ProjectionMatrix <: ProjectiveEntity
end

struct  ProjectionMatrices{T₁ <: ProjectiveEntity, T₂ <: FactorisationAlgorithm, T₃ <: TotalViews } <: ProjectiveEntity
    src::T₁
    method::T₂
    views::T₃
end

struct  LatentVariables{T₁ <: ProjectiveEntity} <: ProjectiveEntity
    src::T₁
end

struct  Chojnacki <: FactorisationAlgorithm
end

mutable struct  HomogeneousCoordinates <: ProjectiveEntity
end

struct CartesianSystem <: PlanarCoordinateSystem
    𝐞₁::Vec{2,Float64}
    𝐞₂::Vec{2,Float64}
end
# Convention I:
# CartesianSystem() = CartesianSystem(Vec(1.0, 0.0), Vec(0.0, -1.0))
CartesianSystem() = CartesianSystem(Vec(-1.0, 0.0), Vec(0.0, 1.0))

struct RasterSystem <: PlanarCoordinateSystem
    𝐞₁::Vec{2,Float64}
    𝐞₂::Vec{2,Float64}
end
#TODO Will have to change this to accomodate different conventions.
# Convention I:  RasterSystem() = RasterSystem(Vec(1.0, 0.0), Vec(0.0, 1.0))
RasterSystem() = RasterSystem(Vec(-1.0, 0.0), Vec(0.0, -1.0))



struct OpticalSystem <: PlanarCoordinateSystem
    𝐞₁::Vec{2,Float64}
    𝐞₂::Vec{2,Float64}
end
# Convention I: = OpticalSystem(Vec(1.0, 0.0), Vec(0.0, 1.0))
OpticalSystem() = OpticalSystem(Vec(-1.0, 0.0), Vec(0.0, -1.0))

mutable struct  Pinhole <: CameraModel
    image_width::Int
    image_height::Int
    focal_length::Float64
    # Center of projection.
    𝐜::Point{3,Float64}
    # Basis vectors that characterise the pose of the camera
    𝐞₁::Vec{3,Float64}
    𝐞₂::Vec{3,Float64}
    𝐞₃::Vec{3,Float64}
    # Origin of the picture plane (the image).
    𝐨::Point{2,Float64}
    # Basis vectors that characterise the coordinate system of the
    # picture plane (the image).
    𝐞₁′::Vec{2,Float64}
    𝐞₂′::Vec{2,Float64}
end

struct  CanonicalLens <: CameraModel
end

struct DirectLinearTransform <: EstimationAlgorithm
end

struct ManualEstimation{A<:AbstractArray} <: EstimationAlgorithm
    𝚹::A
end

struct BundleAdjustment{EA <: EstimationAlgorithm} <: EstimationAlgorithm
    seed::EA
    max_iter::Int64
    toleranceθ::Float64
end

struct FundamentalNumericalScheme{EA <: EstimationAlgorithm} <: EstimationAlgorithm
    seed::EA
    max_iter::Int64
    toleranceθ::Float64
end

mutable struct CanonicalApproximation <: CovarianceEstimationScheme
end

mutable struct HessianApproximation <: CovarianceEstimationScheme
end

mutable struct Taubin <: EstimationAlgorithm
end

mutable struct ApproximateMaximumLikelihood <: CostFunction
end
const AML = ApproximateMaximumLikelihood

struct ReprojectionError <: CostFunction
end

mutable struct CanonicalToHartley <: CoordinateSystemTransformation
end

mutable struct HartleyToCanonical <: CoordinateSystemTransformation
end

mutable struct CovarianceMatrices
end


mutable struct GaussianNoise <: NoiseModel

end

mutable struct TwoViews <: TotalViews
end

abstract type AbstractExperiment end
abstract type AbstractExperimentTrial end
abstract type AbstractExperimentCondition end
abstract type AbstractParticipant end
abstract type AbstractTrialResult end
abstract type AbstractConditionResult end
abstract type AbstractParticipantResult end
abstract type AbstractExperimentResult end


struct Experiment{ActualProjectiveEntity <: ProjectiveEntity, ActualParticipant <: AbstractParticipant, ActualCostFunction <: CostFunction} <: AbstractExperiment
    description::String
    task::ActualProjectiveEntity
    participants::Vector{ActualParticipant}
    cost_function::ActualCostFunction
end

struct ExperimentTrial{ActualArray₁ <: AbstractArray, ActualArray₂ <: AbstractArray , ActualArray₃ <: AbstractArray, ActualArray₄ <: AbstractArray} <: AbstractExperimentTrial
    pure_training_data::Tuple{ActualArray₁, Vararg{ActualArray₁}}
    perturbed_training_data::Tuple{ActualArray₂, Vararg{ActualArray₂}}
    pure_testing_data::Tuple{ActualArray₃, Vararg{ActualArray₃}}
    correct_parameters::ActualArray₄
end

struct ExperimentCondition{ActualExperimentTrial <: AbstractExperimentTrial} <: AbstractExperimentCondition
    description::String
    trials::Vector{ActualExperimentTrial}
end

struct Participant{ActualEstimationAlgorithm <: EstimationAlgorithm, ActualExperimentCondition <: AbstractExperimentCondition} <: AbstractParticipant
    algorithm::ActualEstimationAlgorithm
    conditions::Dict{String, ActualExperimentCondition}
end

struct TrialResult{ActualCostFunction <: CostFunction} <: AbstractTrialResult
    cost_function::ActualCostFunction
    total_points::Int
    val::Vector{Float64}
end

struct ConditionResult{ActualTrialResult <: AbstractTrialResult} <: AbstractConditionResult
    description::String
    outcomes::Vector{ActualTrialResult}
end

struct ParticipantResult{ActualParticipant <: AbstractParticipant, ActualConditionResult <: AbstractConditionResult} <: AbstractParticipantResult
    participant::ActualParticipant
    conditions::Dict{String, ActualConditionResult}
    rmse::Float64
end

struct ExperimentResult{ActualProjectiveEntity <: ProjectiveEntity, ActualParticipantResult <: AbstractParticipantResult, ActualCostFunction <: CostFunction} <: AbstractExperimentResult
    description::String
    task::ActualProjectiveEntity
    results::Vector{ActualParticipantResult}
    cost_function::ActualCostFunction
end

struct PlanarScene
    plane_count::Int
    max_point_count::Int
    roi_range::Tuple{UnitRange{Int},UnitRange{Int}}
end

#Tuple{MeasurementProperties, Vararg{MeasurementProperties}}
