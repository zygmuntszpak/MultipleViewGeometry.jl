abstract type AbstractCoordinateSystem end
abstract type AbstractPlanarCoordinateSystem <: AbstractCoordinateSystem end

Base.@kwdef struct PlanarCartesianSystem{T₁ <: AbstractVector, T₂ <: AbstractVector} <: AbstractPlanarCoordinateSystem
    𝐨::T₁ = Vec(0.0, 0.0)
    𝐞₁::T₂ = Vec(-1.0, 0.0)
    𝐞₂::T₂ = Vec(0.0, 1.0)
end

# TODO Document left-handed versus right-handed convention
Base.@kwdef struct CartesianSystem{T₁ <: AbstractVector, T₂ <: AbstractVector} <: AbstractCoordinateSystem
    𝐨::T₁ = Vec(0.0, 0.0, 0.0)
    𝐞₁::T₂ = Vec(1.0, 0.0, 0.0)
    𝐞₂::T₂ = Vec(0.0, 1.0, 0.0)
    𝐞₃::T₂ = Vec(0.0, 0.0, 1.0)
end

Base.@kwdef struct RasterSystem{T₁ <: AbstractVector, T₂ <: AbstractVector} <: AbstractPlanarCoordinateSystem
    𝐨::T₁ = Vec(0.0, 0.0)
    𝐞₁::T₂ = Vec(-1.0, 0.0)
    𝐞₂::T₂ = Vec(0.0, -1.0)
end

Base.@kwdef struct OpticalSystem{T₁ <: AbstractVector, T₂ <: AbstractVector} <: AbstractPlanarCoordinateSystem
    𝐨::T₁ = Vec(0.0, 0.0)
    𝐞₁::T₂ = Vec(-1.0, 0.0)
    𝐞₂::T₂ = Vec(0.0, -1.0)
end

function origin(coordinate_system::AbstractCoordinateSystem)
    @unpack 𝐨 = coordinate_system
    return 𝐨
end

function basis_vectors(coordinate_system::AbstractCoordinateSystem)
    @unpack 𝐞₁, 𝐞₂, 𝐞₃ = coordinate_system
    return 𝐞₁, 𝐞₂, 𝐞₃
end

function basis_vectors(coordinate_system::AbstractPlanarCoordinateSystem)
    @unpack 𝐞₁, 𝐞₂ = coordinate_system
    return 𝐞₁, 𝐞₂
end

# function origin(param::AbstractCoordinateSystem)
#     param.𝐨
# end
#
# function get_e₁(param::AbstractCoordinateSystem)
#     param.𝐞₁
# end
#
# function get_e₂(param::AbstractCoordinateSystem)
#     param.𝐞₂
# end
#
# function get_e₃(param::AbstractCoordinateSystem)
#     param.𝐞₃
# end
