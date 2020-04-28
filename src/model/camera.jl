abstract type AbstractCamera end
abstract type AbstractCameraModel end
abstract type AbstractIntrinsicParameters end
abstract type AbstractExtrinsicParameters end


Base.@kwdef struct IntrinsicParameters <: AbstractIntrinsicParameters
    focal_length::Float64 = 50
    width::Int = 1000
    height::Int = 1000
    # Diagonal distortion of the image plane (which is usually negligible or zero).
    skewedness::Float64 = 0.0
    # Possibily different sensor scale in the x and y direction, respectively.
    scale_x::Float64 = 1.0
    scale_y::Float64 = 1.0
    # The principal point offset with respect to the Optical Axis coordinate system.
    principal_point::Point{2,Float64} = Point(0.0, 0.0)
    # Origin and basis vectors that characterise the coordinate system of the
    # picture plane (the image).
    coordinate_system = OpticalSystem()
end

function matrix(intrinsics::IntrinsicParameters)
    @unpack focal_length, skewedness, scale_x, scale_y, principal_point = intrinsics
    f = focal_length
    sx = scale_x
    sy = scale_y
    s = skewedness
    𝐩 = principal_point
    𝐊 = @SMatrix [f*sx f*s 𝐩[1]; 0 f*sy 𝐩[2]; 0 0 1]
    return 𝐊
end

Base.@kwdef struct ExtrinsicParameters <: AbstractExtrinsicParameters
    # Origin and basis vectors that characterise the pose of the camera.
    coordinate_system = CartesianSystem(Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
end

function origin(param::ExtrinsicParameters)
    @unpack coordinate_system = param
    return origin(coordinate_system)
end

function basis_vectors(param::ExtrinsicParameters)
    @unpack coordinate_system = param
    return basis_vectors(coordinate_system)
end

Base.@kwdef struct  Pinhole{T₁ <: AbstractIntrinsicParameters, T₂ <: AbstractExtrinsicParameters} <: AbstractCameraModel
    intrinsics::T₁ = IntrinsicParameters()
    extrinsics::T₂ = ExtrinsicParameters()
end

function intrinsics(model::AbstractCameraModel)
    @unpack intrinsics = model
    return intrinsics
end

function extrinsics(model::AbstractCameraModel)
    @unpack extrinsics = model
    return extrinsics
end

Base.@kwdef struct Camera{T₁ <: AbstractCameraModel, T₂ <: AbstractImage} <: AbstractCamera
    model::T₁ = Pinhole()
    image_type::T₂ = AnalogueImage()
end

function model(camera::AbstractCamera)
    @unpack model = camera
    return model
end

function image_type(camera::AbstractCamera)
    @unpack image_type = camera
    return image_type
end

function rotate(camera::Camera,  𝐑::AbstractArray)
    @unpack model = camera
    @unpack extrinsics = model
    @unpack coordinate_system = extrinsics
    @unpack 𝐨, 𝐞₁, 𝐞₂, 𝐞₃ =  coordinate_system
    𝐞₁′ = 𝐑*𝐞₁
    𝐞₂′ = 𝐑*𝐞₂
    𝐞₃′ = 𝐑*𝐞₃
    return @set camera.model.extrinsics.coordinate_system = CartesianSystem(𝐨, 𝐞₁′,𝐞₂′,𝐞₃′)
end

function translate(camera::Camera, 𝐭::AbstractArray)
    @unpack model = camera
    @unpack extrinsics = model
    @unpack coordinate_system = extrinsics
    @unpack 𝐨, 𝐞₁, 𝐞₂, 𝐞₃ =  coordinate_system
    return @set camera.model.extrinsics.coordinate_system = CartesianSystem(𝐨 + 𝐭, 𝐞₁, 𝐞₂, 𝐞₃)
end

function relocate(camera::Camera, 𝐑::AbstractArray, 𝐭::AbstractArray)
    @unpack model = camera
    @unpack extrinsics = model
    @unpack coordinate_system = extrinsics
    @unpack 𝐨, 𝐞₁, 𝐞₂, 𝐞₃ =  coordinate_system
    𝐞₁′ = 𝐑*𝐞₁
    𝐞₂′ = 𝐑*𝐞₂
    𝐞₃′ = 𝐑*𝐞₃
    return @set camera.model.extrinsics.coordinate_system = CartesianSystem(𝐨 + 𝐭, 𝐞₁′, 𝐞₂′, 𝐞₃′)
end
