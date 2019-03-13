"""
    ascertain_pose(camera::CameraModel, 𝐞₁, 𝐞₂, 𝐞₃)

Determines the rotation and translation of the camera with respect to the
origin of a world coordinate system with axes {𝐞₁, 𝐞₂, 𝐞₃}.
"""
function ascertain_pose(camera::CameraModel, 𝐞₁, 𝐞₂, 𝐞₃)
    𝐭 = camera.𝐜
    𝐑 = inv(hcat(𝐞₁, 𝐞₂, 𝐞₃)) * hcat(camera.𝐞₁, camera.𝐞₂, camera.𝐞₃)
    𝐑, 𝐭
end

function obtain_intrinsics(camera::CameraModel, system::RasterSystem)
    # Convention I: 𝐭 = Point(-camera.image_width / 2, camera.image_height / 2)
    𝐭 = Point(-camera.image_width / 2, -camera.image_height / 2)
    𝐑 = inv(hcat(system.𝐞₁, system.𝐞₂)) * hcat(camera.𝐞₁′, camera.𝐞₂′)
    f = camera.focal_length
    𝐊 = SMatrix{3,3,Float64,9}(f, 0.0, 0.0, 0.0, f, 0.0, 0.0, 0.0 , 1)
    𝐊′ =vcat(hcat(𝐑', -𝐑'*𝐭), SMatrix{1,3,Float64}(0,0,1) )*𝐊
end

function obtain_intrinsics(camera::CameraModel, system::CartesianSystem)
    # Convention I: 𝐭 = Point(-camera.image_width / 2, camera.image_height / 2)
    𝐭 = Point(-camera.image_width / 2, camera.image_height / 2)
    𝐑 = inv(hcat(system.𝐞₁, system.𝐞₂)) * hcat(camera.𝐞₁′, camera.𝐞₂′)
    #𝐑 = vcat(hcat(𝐑₂₂, 0), SVector(0,0,0))
    #𝐑 =vcat(hcat(𝐑₂₂, SVector(0,0)), SMatrix{1,3,Float64}(0,0,1))
    f = camera.focal_length
    𝐊 = SMatrix{3,3,Float64,9}(f, 0.0, 0.0, 0.0, f, 0.0, 0.0, 0.0 , 1)
    𝐊′ =vcat(hcat(𝐑', -𝐑'*𝐭), SMatrix{1,3,Float64}(0,0,1))*𝐊
end

function obtain_intrinsics(camera::CameraModel, system::OpticalSystem)
    𝐭 = Point(0, 0)
    𝐑 = inv(hcat(system.𝐞₁, system.𝐞₂)) * hcat(camera.𝐞₁′, camera.𝐞₂′)
    @show 𝐑
    f = camera.focal_length
    𝐊 = SMatrix{3,3,Float64,9}(f, 0.0, 0.0, 0.0, f, 0.0, 0.0, 0.0 , 1)
    𝐊′ =vcat(hcat(𝐑', -𝐑'*𝐭), SMatrix{1,3,Float64}(0,0,1))*𝐊
end
