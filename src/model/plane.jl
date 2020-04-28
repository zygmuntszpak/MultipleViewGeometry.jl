abstract type AbstractPlane end

Base.@kwdef struct Plane{T₁ <: AbstractVector, T₂ <: Real} <: AbstractPlane
    normal::T₁
    distance::T₂
end

function normal(plane::Plane)
    @unpack normal = plane
    return normal
end

function distance(plane::Plane)
     @unpack distance = plane
    return distance
end

Base.@kwdef struct PlaneSegment{T₁ <: Plane, T₂ <: AbstractVector} <: AbstractPlane
    plane::T₁
    segment::T₂
end

function plane(plane_segment::PlaneSegment)
    @unpack plane = plane_segment
    return plane
end

# function set_plane!(plane_segment::PlaneSegment, plane::Plane)
#     plane_segment.plane = plane
# end

function normal(plane_segment::PlaneSegment)
    return normal(plane(plane_segment))
end

function distance(plane_segment::PlaneSegment)
    return distance(plane(plane_segment))
end

function segment(plane_segment::PlaneSegment)
    @unpack segment = plane_segment
    return segment
end

function on_plane(𝐗::AbstractVector, plane::Plane; tol::Number = 1e-10)
    𝐧 = normal(plane)
    d = distance(plane)
    return abs(dot(𝐗,𝐧) - d) < tol ? true : false
end
