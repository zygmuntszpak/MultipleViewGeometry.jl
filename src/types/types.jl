

struct HomogeneousPoint{T <: AbstractFloat,N}
    coords::NTuple{N, T}
end

abstract type ProjectiveEntity end

abstract type CameraModel end

abstract type EstimationAlgorithm end

abstract type CostFunction end

abstract type CoordinateSystemTransformation end

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
