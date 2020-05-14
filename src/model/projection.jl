struct Projection{T₁ <: AbstractMatrix} <: ProjectiveEntity
    𝐏::T₁
end

function matrix(entity::Projection)
    entity.𝐏
end

Projection(camera::AbstractCamera) = Projection(construct_projection(camera, CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))))
Projection(camera::AbstractCamera, world_coordinate_system::AbstractCoordinateSystem) = Projection(construct_projection(camera, world_coordinate_system))
function project(P::Projection, 𝒳::Vector{<: AbstractVector})
    𝐏 = matrix(P)
    ℳ = map(𝒳) do 𝐗
        𝐦 = hom⁻¹(𝐏 * hom(𝐗))
    end
    return ℳ
end

function back_project(camera::AbstractCamera, ℳ::Vector{<: AbstractVector})
    @unpack model, image_type = camera
    @unpack extrinsics, intrinsics = model
    @unpack 𝐨, 𝐞₁, 𝐞₂, 𝐞₃ = extrinsics
    @unpack coordinate_system = image_type
    @unpack focal_length = intrinsics
    ℒ = map(ℳ) do 𝐦
        𝐩 = 𝐨 + 𝐦[1] *𝐞₁ + 𝐦[2] *𝐞₂ + focal_length*𝐞₃
        L = Line3D(𝐨, 𝐩)
    end
    return ℒ
end

function construct_projection(camera::AbstractCamera, world_reference_system::AbstractCoordinateSystem)
    @unpack model, image_type = camera
    @unpack coordinate_system = image_type
    construct_projection(model, world_reference_system, coordinate_system)
end

# World reference system is always identity
function construct_projection(model::AbstractCameraModel, world_reference_system::AbstractCoordinateSystem, image_system::AbstractPlanarCoordinateSystem)
    @unpack intrinsics, extrinsics = model
    𝐊 = matrix(intrinsics, image_system)
    𝐄 = matrix(extrinsics, world_reference_system)
    𝐏 = 𝐊 * 𝐄
    𝐏 = 𝐏 / norm(𝐏) # TODO optionally remove this normalization
    return 𝐏
end

function matrix(intrinsics::IntrinsicParameters, image_system::AbstractPlanarCoordinateSystem)
    @unpack focal_length, coordinate_system  = intrinsics
    @unpack skewedness, scale_x, scale_y, principal_point = intrinsics
    𝐞₁, 𝐞₂ = basis_vectors(image_system)
    𝐞₁′, 𝐞₂′ = basis_vectors(coordinate_system)
    f = focal_length
    sx = scale_x
    sy = scale_y
    s = skewedness
    𝐩 = principal_point
    # TODO Fix this so that we don't assume that the principal point is at position (0,0)
    𝐭 = determine_translation(intrinsics, image_system)
    @warn "Need to ensure this `matrix` function takes into account non-zero principal point"
    𝐑 = inv(hcat(𝐞₁, 𝐞₂)) * hcat(𝐞₁′ , 𝐞₂′)
    𝐊 = @SMatrix [f*sx f*s 𝐩[1]; 0 f*sy 𝐩[2]; 0 0 1]
    𝐊′ = vcat(hcat(𝐑', -𝐑'*𝐭), SMatrix{1,3,Float64}(0,0,1)) * 𝐊
    return 𝐊′
end

function matrix(extrinsics::ExtrinsicParameters, world_reference_system::CartesianSystem = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0)))
    𝐑, 𝐭 = ascertain_pose(extrinsics, world_reference_system)
    𝐄 = [𝐑' -𝐑'*𝐭] # TODO convert to static array
    return 𝐄
end

function ascertain_pose(camera::AbstractCamera, reference_system::CartesianSystem)
    @unpack model = camera
    @unpack extrinsics = model
    return ascertain_pose(extrinsics, reference_system)
end

function ascertain_pose(extrinsics::ExtrinsicParameters, reference_system::CartesianSystem)
    @unpack coordinate_system = extrinsics
    𝐞₁, 𝐞₂, 𝐞₃ = basis_vectors(reference_system)
    𝐞₁′, 𝐞₂′, 𝐞₃′ = basis_vectors(coordinate_system)
    𝐭 = origin(coordinate_system) - origin(reference_system)
    𝐑 = inv(hcat(𝐞₁, 𝐞₂, 𝐞₃)) * hcat(𝐞₁′, 𝐞₂′, 𝐞₃′)
    return 𝐑, 𝐭
end

# TODO Incorporate information about the origin of the image coordinate system
function determine_translation(intrinsics::IntrinsicParameters, system::PlanarCartesianSystem)
    @unpack width, height = intrinsics
    𝐭 = Point(-width / 2, height / 2)
    return 𝐭
end

# TODO Incorporate information about the origin of the image coordinate system
function determine_translation(intrinsics::IntrinsicParameters, system::OpticalSystem)
    𝐭 = Point(0, 0)
    return 𝐭
end

# TODO Incorporate information about the origin of the image coordinate system
function determine_translation(intrinsics::IntrinsicParameters, system::RasterSystem)
    @unpack width, height = intrinsics
    𝐭 = Point(-width / 2, -height / 2)
    return 𝐭
end
