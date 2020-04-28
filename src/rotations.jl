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

function rodrigues2matrix(vˣ::Real,vʸ::Real,vᶻ::Real)
      𝐯 = [vˣ, vʸ, vᶻ]
      θ = norm(𝐯)
      𝐯 = θ == 0 ? 𝐯 : 𝐯/θ
      𝐈 = Matrix(1.0I, 3, 3)
      𝐖 = vec2antisym(𝐯)
      𝐑 = 𝐈 + 𝐖 * sin(θ) + 𝐖^2 * (1-cos(θ))
end
