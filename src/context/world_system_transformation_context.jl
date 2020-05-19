abstract type AbstractCoordinateTransformationContext <: AbstractContext end

struct WorldSystemTransformation{T <: CoordinateTransformation} <: AbstractCoordinateTransformationContext
     coordinate_transformation::T
end

function (context::WorldSystemTransformation)(camera::AbstractCamera)
    𝐑 = rotation(context.coordinate_transformation)
    𝐭 = translation(context.coordinate_transformation)
    @unpack coordinate_system = extrinsics(camera)
    𝐞₁, 𝐞₂, 𝐞₃ = basis_vectors(coordinate_system)
    𝐨 = origin(coordinate_system)
    𝐞₁′ = 𝐑' * 𝐞₁
    𝐞₂′ = 𝐑' * 𝐞₂
    𝐞₃′ = 𝐑' * 𝐞₃
    𝐨′ =  𝐑' * (𝐨 - 𝐭)
    return @set camera.model.extrinsics.coordinate_system = CartesianSystem(𝐨′, 𝐞₁′, 𝐞₂′, 𝐞₃′)
end

function (context::WorldSystemTransformation)(world::AbstractWorld)
    @unpack coordinate_transformation = context
    @unpack points, planes = world
    𝐑 = rotation(coordinate_transformation)
    𝐭 = translation(coordinate_transformation)
    points′ = transform_3D_points(𝐑, 𝐭, points)
    planes′ = transform_planes(𝐑, 𝐭, planes)
    world = @set world.points = points′
    world = @set world.planes = planes′
    return world
end

# TODO extend functor for planes and points.

function transform_3D_points(𝐑::AbstractMatrix, 𝐭::AbstractVector, points::AbstractVector)
    map(points) do 𝐗
        𝐑' * (𝐗 - 𝐭)
    end
end

function transform_planes(𝐑::AbstractMatrix, 𝐭::AbstractVector, planes::Vector{<: Union{Plane, PlaneSegment}})
    [transform_plane(𝐑, 𝐭, planes[k]) for k = 1:length(planes)]
end

function transform_plane(𝐑::AbstractMatrix, 𝐭::AbstractVector, plane::T) where T <: Union{Plane, PlaneSegment}
    𝐧 = normal(plane)
    d = distance(plane)
    𝐚 = construct_point_on_plane(𝐧, d)

    𝐧′ = 𝐑' * 𝐧
    𝐚′ = 𝐑' * (𝐚 - 𝐭)
    d′ = dot(𝐚′, 𝐧′)

    # Ensure that our plane representation always follows the "outward normal" convention
    if d′ < 0
        return T(Vec3(-𝐧′...), -d′)
    else
        return T(Vec3(𝐧′...), d′)
    end
end

function transform_planes(𝐑::AbstractMatrix, 𝐭::AbstractVector, planes::Vector{<: EuclideanPlane3D})
    [transform_plane(𝐑, 𝐭, planes[k]) for k = 1:length(planes)]
end

function transform_plane(𝐑::AbstractMatrix, 𝐭::AbstractVector, plane::T) where T <: Union{EuclideanPlane3D}
    @unpack  coordinate_system = plane
    𝐞₁, 𝐞₂, 𝐞₃ = basis_vectors(coordinate_system)
    𝐨 = origin(coordinate_system)
    𝐞₁′ = 𝐑' * 𝐞₁
    𝐞₂′ = 𝐑' * 𝐞₂
    𝐞₃′ = 𝐑' * 𝐞₃
    𝐨′ =  𝐑' * (𝐨 - 𝐭)
    plane′ = EuclideanPlane3D(CartesianSystem(𝐨′, 𝐞₁′, 𝐞₂′, 𝐞₃′))
    d′ = distance(plane′)
    if d′ < 0
        return EuclideanPlane3D(CartesianSystem(-𝐨′, 𝐞₁′, 𝐞₂′, -𝐞₃′))
    else
        return plane′
    end
end

function construct_point_on_plane(𝐧::AbstractVector, d::Number)
    if 𝐧[1] != 0
        a = d / 𝐧[1]
        𝐚 = [a, 0.0, 0.0]
        return 𝐚
    elseif 𝐧[2] != 0
        a = d / 𝐧[2]
        𝐚 = [0.0, a, 0.0]
        return 𝐚
    else
        a = d / 𝐧[3]
        𝐚 = [0.0, 0.0, a]
        return 𝐚
    end
end
