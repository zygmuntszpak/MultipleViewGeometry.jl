abstract type AbstractImage end
abstract type AbstractAnalogueImage <: AbstractImage end
abstract type AbstractDigitalImage <: AbstractImage end

Base.@kwdef struct AnalogueImage{T₁ <: AbstractPlanarCoordinateSystem} <: AbstractAnalogueImage
    coordinate_system::T₁ = PlanarCartesianSystem(𝐞₁ = Vec(-1.0, 0.0), 𝐞₂ = Vec(0.0, -1.0))
end

function coordinate_system(image::AbstractAnalogueImage)
    @unpack coordinate_system = image
    return coordinate_system
end

# function get_coordinate_system(image::AnalogueImage)
#     image.coordinate_system
# end

# function get_data(image::AnalogueImage)
#     image.data
# end
