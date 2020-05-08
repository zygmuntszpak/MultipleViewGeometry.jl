"""julia
    rotx(θ::Real)

A matrix representing the rotation about the X-axis by an angle θ.
"""
function rotx(θ::Real)
cosθ = cos(θ);
sinθ = sin(θ);
𝐑ˣ = [1.0  0.0      0.0  ;
      0.   cosθ    -sinθ ;
      0    sinθ     cosθ ]
end

"""juliaa
    roty(θ::Real)

A matrix representing rotation about the Y-axis by an angle θ.
"""
function roty(θ::Real)
cosθ = cos(θ);
sinθ = sin(θ);
𝐑ʸ = [cosθ     0.0    sinθ ;
      0.0      1.0    0.0  ;
      -sinθ    0.0    cosθ]
end

"""julia
    rotz(θ::Real)

A matrix representing the rotation about the Z-axis by an angle θ.
"""
function rotz(θ::Real)
cosθ = cos(θ);
sinθ = sin(θ);
𝐑ᶻ = [cosθ     -sinθ    0.0 ;
      sinθ      cosθ    0.0 ;
      0.0       0.0     1.0 ]
end

"""julia
    rotz(θ::Real)

A matrix representing the rotation about the X-axis, Y-axis and Z-axis by the angles θˣ, θʸ, and θᶻ respectively.
"""
function rotxyz(θˣ::Real,θʸ::Real,θᶻ::Real)
      rotx(θˣ)*roty(θʸ)*rotz(θᶻ)
end

# function rodrigues2matrix(vˣ::Real,vʸ::Real,vᶻ::Real)
#       𝐯 = [vˣ, vʸ, vᶻ]
#       θ = norm(𝐯)
#       𝐯 = θ == 0 ? 𝐯 : 𝐯/θ
#       𝐈 = Matrix(1.0I, 3, 3)
#       𝐖 = vec2antisym(𝐯)
#       𝐑 = 𝐈 + 𝐖 * sin(θ) + 𝐖^2 * (1-cos(θ))
# end

# """
#     matrix2rodrigues(𝐑::AbstratArray)
#
#     Conversion between a rotation matrix and the associated Rodrigues rotation vector.
#
#
# # Reference
# [1] Burger, Wilhelm. (2016). Zhang's Camera Calibration Algorithm: In-Depth Tutorial and Implementation.
# [2] C. Tomasi. Vector representation of rotations. Computer Science 527
# Course Notes, Duke University, https://www.cs.duke.edu/courses/fall13/
# compsci527/notes/rodrigues.pdf, 2013.
# """
# function matrix2rodrigues(𝐑::AbstractArray)
#       𝐑₃₂ =
#       𝐩 = 0.5 * SVector(𝐑[3,2] - 𝐑[2,3], 𝐑[1,3] - 𝐑[3,1], 𝐑[2,1] - 𝐑[1,2])
#       c = 0.5 * (trace(𝐑) - 1.0)
# end
