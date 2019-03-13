function construct( e::HomographyMatrix,
                   𝐊₁::AbstractArray,
                   𝐑₁::AbstractArray,
                   𝐭₁::AbstractArray,
                   𝐊₂::AbstractArray,
                   𝐑₂::AbstractArray,
                   𝐭₂::AbstractArray,
                   𝐧::AbstractArray,
                   d::Real)

    if size(𝐊₁) != (3,3) || size(𝐊₂) != (3,3) ||
       size(𝐑₁) != (3,3) || size(𝐑₂) != (3,3)
        throw(ArgumentError("Expect 3 x 3 calibration and rotation matrices."))
    end
    if length(𝐭₁) != 3 || length(𝐭₂) != 3
        throw(ArgumentError("Expect length-3 translation vectors."))
    end
    if length(𝐧) != 3
        throw(ArgumentError("Expect length-3 normal vector."))
    end
    # TODO Check that camera center does not lie on the plane. 
    𝐈 = SMatrix{3,3}(1.0I)
    𝐇 = 𝐊₂*𝐑₂*(𝐈 - (𝐧'*𝐭₁ + d)^-1 * (𝐭₁ - 𝐭₂) * 𝐧')/𝐑₁/𝐊₁
    SMatrix{3,3,Float64,3*3}(𝐇)
end
