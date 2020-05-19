Base.@kwdef struct EuclideanPlane3D <: AbstractPlane
    coordinate_system::CartesianSystem = CartesianSystem()
end

# TODO include constructor for
#EuclideanPlane3D(origin::AbstractVector, normal::AbstractVector) = ...


function normal(plane::EuclideanPlane3D)
    @unpack coordinate_system = plane
    @unpack 𝐞₃ = coordinate_system
    return 𝐞₃
end

function origin(plane::EuclideanPlane3D)
    @unpack coordinate_system = plane
    return origin(coordinate_system)
end


function distance(plane::EuclideanPlane3D)
    dot(normal(plane), origin(plane))
end

function on_plane(𝐗::AbstractVector, plane::EuclideanPlane3D; tol::Number = 1e-10)
    𝐧 = normal(plane)
    d = distance(plane)
    abs(dot(𝐗,𝐧) - d) < tol ? true : false
end

Base.@kwdef struct Line3D{T <: AbstractVector} <: AbstractPlane
    𝐩₁::T = Vec(0.0, 0.0, 0.0)
    𝐩₂::T = Vec(0.0, 0.0, 1.0)
end

function on_line(𝐩::AbstractVector, 𝓁::Line3D; tol = 1e-10)
    𝐯₁ =  𝓁.𝐩₁ - 𝐩
    𝐯₂ =  𝓁.𝐩₂ - 𝐩
    n = norm(cross(𝐯₁, 𝐯₂))
    return  n <= tol ? true : false
end
