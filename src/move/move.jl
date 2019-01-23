"""
    translate!(camera::CameraModel, 𝐭::AbstractArray)

Translate the camera optical center by a vector 𝐭.
"""
function translate!(camera::CameraModel, 𝐭::AbstractArray)
    camera.𝐜 = camera.𝐜 + 𝐭
    camera
end

"""
    relocate!(camera::CameraModel, 𝐭::AbstractArray)

Rotates a camera around its optical center by a rotation matrix 𝐑 and
then translates the optical center by a vector 𝐭.
"""
function relocate!(camera::CameraModel, 𝐑::AbstractArray, 𝐭::AbstractArray)
    camera.𝐞₁ = 𝐑*camera.𝐞₁
    camera.𝐞₂ = 𝐑*camera.𝐞₂
    camera.𝐞₃ = 𝐑*camera.𝐞₃
    camera.𝐜 = camera.𝐜 + 𝐭
    camera
end
