abstract type AbstractPose end
abstract type AbstractCoordinateTransformation end

struct RelativePose{T₁ <: AbstractMatrix, T₂ <: AbstractVector} <: AbstractPose
    rotation::T₁
    translation::T₂
end

function rotation(pose::AbstractPose)
    @unpack rotation = pose
    return rotation
end

function translation(pose::AbstractPose)
    @unpack translation = pose
    return translation
end

function matrix(pose::AbstractPose)
    @unpack rotation, translation = pose
    𝐄 = hcat(rotation, translation)
    return 𝐄
end

RelativePose(src_camera::AbstractCamera, tgt_camera::AbstractCamera) = RelativePose(construct_relative_pose(src_camera, tgt_camera)...)
RelativePose(coordinate_system₁::CartesianSystem, coordinate_system₂::CartesianSystem) = RelativePose(construct_relative_pose(coordinate_system₁, coordinate_system₂)...)

# tgt_camera with respect to src_camera
function construct_relative_pose(src_camera::AbstractCamera, tgt_camera::AbstractCamera)
    src_extrinsics = extrinsics(src_camera)
    tgt_extrinsics = extrinsics(tgt_camera)
    src_camera_system = coordinate_system(src_extrinsics)
    tgt_camera_system = coordinate_system(tgt_extrinsics)
    construct_relative_pose(src_camera_system₁, tgt_camera_system₂)
end

# tgt_system with respect to src_system
function construct_relative_pose(src_coordinate_system::AbstractCoordinateSystem,  tgt_coordinate_system::AbstractCoordinateSystem)
    𝐞₁, 𝐞₂, 𝐞₃ = basis_vectors(src_coordinate_system)
    𝐞₁′, 𝐞₂′, 𝐞₃′ = basis_vectors(tgt_coordinate_system)
    # ? Mistake, change order TODO ?
    𝐭 = origin(tgt_coordinate_system) - origin(src_coordinate_system)
    𝐑 = inv(hcat(𝐞₁, 𝐞₂, 𝐞₃)) * hcat(𝐞₁′, 𝐞₂′, 𝐞₃′)
    𝐑, 𝐭
end

Base.@kwdef struct CoordinateTransformation{T₁ <: AbstractCoordinateSystem, T₂ <: AbstractCoordinateSystem} <: AbstractCoordinateTransformation
    source::T₁ = CartesianSystem()
    target::T₂ = CartesianSystem()
    relative_pose = RelativePose(source, target)
end

function rotation(transformation::CoordinateTransformation)
    @unpack relative_pose = transformation
    return rotation(relative_pose)
end

function translation(transformation::CoordinateTransformation)
    @unpack relative_pose = transformation
    return translation(relative_pose)
end
