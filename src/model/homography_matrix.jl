struct HomographyMatrix{T₁ <: AbstractMatrix} <: ProjectiveEntity
    𝐇::T₁
end

struct HomographyMatrices{T₁ <: AbstractVector{<: HomographyMatrix}} <: ProjectiveEntity
    ℋ::T₁
end

function matrix(entity::HomographyMatrix)
    entity.𝐇
end

function matrices(entity::HomographyMatrices)
    map(x->matrix(x), entity.ℋ)
end

HomographyMatrices(camera₁::AbstractCamera, camera₂::AbstractCamera, planes::AbstractVector{<:Union{Plane, EuclideanPlane3D}}) = HomographyMatrices(construct_homography_matrices(camera₁, camera₂, planes,  CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))))
HomographyMatrix(camera₁::AbstractCamera, camera₂::AbstractCamera, plane::Union{Plane, EuclideanPlane3D}) = HomographyMatrix(construct_homography_matrix(camera₁, camera₂, plane, CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))))

function construct_homography_matrices(camera₁::AbstractCamera, camera₂::AbstractCamera, planes::AbstractVector{<:Union{Plane, EuclideanPlane3D}}, reference_system::AbstractCoordinateSystem)
    ℋ = [HomographyMatrix(construct_homography_matrix(camera₁, camera₂, planes[i], reference_system)) for i = 1:length(planes)]
    return ℋ
end

function construct_homography_matrix(camera₁::AbstractCamera, camera₂::AbstractCamera, plane::Union{Plane, EuclideanPlane3D}, reference_system::AbstractCoordinateSystem)
    model₁ = model(camera₁)
    model₂ = model(camera₂)
    image_type₁ = image_type(camera₁)
    image_system₁ = coordinate_system(image_type₁)
    image_type₂ = image_type(camera₂)
    image_system₂ = coordinate_system(image_type₂)
    return construct_homography_matrix(model₁, model₂, plane, reference_system, image_system₁, image_system₂)
end

function construct_homography_matrix(model₁::AbstractCameraModel, model₂::AbstractCameraModel,  plane::Union{Plane, EuclideanPlane3D},  reference_system::AbstractCoordinateSystem, image_system₁::AbstractPlanarCoordinateSystem, image_system₂::AbstractPlanarCoordinateSystem)
    intrinsics₁ = intrinsics(model₁)
    𝐊₁ = matrix(intrinsics₁, image_system₁)
    extrinsics₁ = extrinsics(model₁)
    𝐑₁′, 𝐭₁′ = ascertain_pose(extrinsics₁, reference_system)
    𝐑₁ = transpose(𝐑₁′)
    𝐭₁ = 𝐭₁′
    # Our projection matrix should decompose as [𝐑 -𝐑*𝐭]

    intrinsics₂ = intrinsics(model₂)
    𝐊₂ = matrix(intrinsics₂, image_system₂)
    extrinsics₂ = extrinsics(model₂)
    𝐑₂′, 𝐭₂′ = ascertain_pose(extrinsics₂, reference_system)
    # Our projection matrix should decompose as [𝐑 -𝐑*𝐭]
    𝐑₂ = transpose(𝐑₂′)
    𝐭₂ = 𝐭₂′

    # We assume that the plane is given by the vector 𝛑 =[n -d], where n is the outward
    # normal to the plane and d is the distance from the plane to the origin of the
    # coordinate system.
    𝐧 = normal(plane)
    d = distance(plane)

    𝐀 = 𝐊₂*𝐑₂/𝐑₁/𝐊₁
    𝐛 = 𝐊₂*𝐑₂*(𝐭₁ - 𝐭₂)
    w = d - 𝐧'*𝐭₁
    𝐯 = inv(𝐊₁')*𝐑₁*𝐧

    𝐇 =  w*𝐀 + 𝐛*𝐯'
    return 𝐇
end
