struct FundamentalMatrix{T₁ <: AbstractMatrix} <: ProjectiveEntity
    𝐅::T₁
end

function matrix(entity::FundamentalMatrix)
    entity.𝐅
end


FundamentalMatrix(camera₁::AbstractCamera, camera₂::AbstractCamera) = FundamentalMatrix(construct_fundamental_matrix(camera₁, camera₂, CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))))
FundamentalMatrix(P₁::Projection, P₂::Projection) = FundamentalMatrix(construct_fundamental_matrix(P₁, P₂))

function construct_fundamental_matrix(camera₁::AbstractCamera, camera₂::AbstractCamera,  world_system_transformation::AbstractCoordinateSystem)
    model₁ = model(camera₁)
    model₂ = model(camera₂)
    image_type₁ = image_type(camera₁)
    image_system₁ = coordinate_system(image_type₁)
    image_type₂ = image_type(camera₂)
    image_system₂ = coordinate_system(image_type₂)
    return construct_fundamental_matrix(model₁, model₂,  world_system_transformation, image_system₁, image_system₂)
end

function construct_fundamental_matrix(model₁::AbstractCameraModel, model₂::AbstractCameraModel,  world_system_transformation::AbstractCoordinateSystem, image_system₁::AbstractPlanarCoordinateSystem, image_system₂::AbstractPlanarCoordinateSystem)
    intrinsics₁ = intrinsics(model₁)
    𝐊₁ = matrix(intrinsics₁, image_system₁)
    extrinsics₁ = extrinsics(model₁)
    𝐑₁′, 𝐭₁′ = ascertain_pose(extrinsics₁, world_system_transformation)
    𝐑₁ = transpose(𝐑₁′)
    𝐭₁ = 𝐭₁′
    # Our projection matrix should decompose as [𝐑 -𝐑*𝐭]

    intrinsics₂ = intrinsics(model₂)
    𝐊₂ = matrix(intrinsics₂, image_system₂)
    extrinsics₂ = extrinsics(model₂)
    𝐑₂′, 𝐭₂′ = ascertain_pose(extrinsics₂, world_system_transformation)
    # Our projection matrix should decompose as [𝐑 -𝐑*𝐭]
    𝐑₂ = transpose(𝐑₂′)
    𝐭₂ = 𝐭₂′

    𝐅 = vec2antisym(𝐊₂*𝐑₂*(𝐭₁ - 𝐭₂))*𝐊₂*𝐑₂/𝐑₁/𝐊₁
    return 𝐅
end

function construct_fundamental_matrix(P₁::Projection, P₂::Projection)
    𝐏₁ = matrix(P₁)
    𝐏₂ = matrix(P₂)
    𝐜₁ = SVector{4,Float64}(nullspace(Array(𝐏₁)))
    𝐞₂ = 𝐏₂*𝐜₁
    𝐅 = vec2antisym(𝐞₂)*𝐏₂*pinv(𝐏₁)
    return SMatrix{3,3,Float64,3*3}(𝐅)
end
