# Generates a random point on the plane centered around a point on the plane
# that is closest to the origin.
function generate_planar_points(𝐧::AbstractArray, d::Real, extent::Real, N::Int)
    # Generate vector 𝐰 on a plane through the origin with normal vector 𝐧.
    first(𝐧) == 0 ? 𝐰 = cross(𝐧,SVector(1.0,0.0,0.0)) : 𝐰 = cross(𝐧,SVector(0.0,0.0,1.0))
    points = Array{SVector{3,Float64},1}(undef,N)
    for n = 1:N
        # Rotate 𝐰 randomly around the axis 𝐧.
        θ = rand() * 2*pi
        𝐤 = 𝐧 / norm(𝐧)
        𝐯 = 𝐰*cos(θ) + cross(𝐤,𝐰)*sin(θ) + 𝐤*dot(𝐤,𝐰)*(1-cos(θ))
        # Scale the vector so that it lies in the interval [0, extent)
        𝐯 = (rand() * extent) * 𝐯
        # Translate the vector so that it lies on the plane parametrised by 𝐧 and d.
        𝐯 = 𝐯 + d*(𝐧/norm(𝐧)^2)
        points[n] = 𝐯
    end
    points
end

function generate_planar_points(x_range::AbstractRange, y_range::AbstractRange, z::Real, N::Int)
    𝒳 = [Point3( rand(x_range), rand(y_range), z) for n = 1:N]
end


function crop(rect::HyperRectangle{2, <:Real}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ₁, ℳ₁ʹ = 𝒟
    ℳ₂ = Array{Point{2,Float64}}(undef,0)
    ℳ₂ʹ= Array{Point{2,Float64}}(undef,0)
    N = length(ℳ₁)
    for n = 1:N
        𝐦 = ℳ₁[n]
        𝐦′ = ℳ₁ʹ[n]
        if first(rect.origin) <= first(𝐦) <= first(rect.origin) + first(rect.widths) &&
           last(rect.origin) <= last(𝐦) <= last(rect.origin) + last(rect.widths)
           push!(ℳ₂, 𝐦)
           push!(ℳ₂ʹ, 𝐦′)
       end
    end
    (ℳ₂, ℳ₂ʹ)
end
