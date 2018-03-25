struct HomogeneousPoint{T <: AbstractFloat,N}
    coords::NTuple{N, T}
end

abstract type ProjectiveEntity end

abstract type CameraModel end

type FundamentalMatrix <: ProjectiveEntity
end

type ProjectionMatrix <: ProjectiveEntity
end

type Pinhole <: CameraModel
end

type CanonicalLens <: CameraModel
end
