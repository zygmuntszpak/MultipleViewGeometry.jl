struct HomogeneousPoint{T <: AbstractFloat,N}
    coords::NTuple{N, T}
end

abstract type ProjectiveEntity end

type FundamentalMatrix <: ProjectiveEntity
end
