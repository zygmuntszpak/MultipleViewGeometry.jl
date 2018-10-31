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

abstract type CostFunction end

abstract type CoordinateSystemTransformation end

abstract type CovarianceEstimationScheme end

abstract type NoiseModel end


mutable struct FundamentalMatrix <: ProjectiveEntity
end

mutable struct EssentialMatrix <: ProjectiveEntity
end

mutable struct  ProjectionMatrix <: ProjectiveEntity
end

mutable struct  HomogeneousCoordinates <: ProjectiveEntity
end

mutable struct  Pinhole <: CameraModel
end

mutable struct  CanonicalLens <: CameraModel
end

mutable struct DirectLinearTransform <: EstimationAlgorithm
end

mutable struct BundleAdjustment{A<:AbstractVector} <: EstimationAlgorithm
    𝛉₀::A
    max_iter::Int64
    toleranceθ::Float64
end

mutable struct FundamentalNumericalScheme{A<:AbstractVector} <: EstimationAlgorithm
    𝛉₀::A
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

mutable struct CanonicalToHartley <: CoordinateSystemTransformation
end

mutable struct HartleyToCanonical <: CoordinateSystemTransformation
end

mutable struct CovarianceMatrices
end


mutable struct GaussianNoise <: NoiseModel

end
