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

const Point2DH = MVector{3,Float64}
const Point3DH = MVector{4,Float64}

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

type FundamentalMatrix <: ProjectiveEntity
end

type ProjectionMatrix <: ProjectiveEntity
end

type HomogeneousCoordinates <: ProjectiveEntity
end

type Pinhole <: CameraModel
end

type CanonicalLens <: CameraModel
end

type DirectLinearTransform <: EstimationAlgorithm
end

type FundamentalNumericalScheme <: EstimationAlgorithm
    𝛉₀::Matrix{Float64}
    max_iter::Int8
    toleranceθ::Float64
end

type CanonicalApproximation <: CovarianceEstimationScheme
end

type HessianApproximation <: CovarianceEstimationScheme
end



type Taubin <: EstimationAlgorithm
end

type ApproximateMaximumLikelihood <: CostFunction
end
const AML = ApproximateMaximumLikelihood

type CanonicalToHartley <: CoordinateSystemTransformation
end

type HartleyToCanonical <: CoordinateSystemTransformation
end

type CovarianceMatrices
end
