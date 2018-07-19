function construct( e::EssentialMatrix, 𝐅::AbstractArray,  𝐊₁::AbstractArray, 𝐊₂::AbstractArray)
    if (size(𝐊₁) != (3,3)) || (size(𝐊₂) != (3,3))
        throw(ArgumentError("Expect 3 x 3 calibration matrices."))
    end
    if (size(𝐅) != (3,3))
        throw(ArgumentError("Expect 3 x 3 fundamental matrix."))
    end
    # Equation 9.12 Chapter 9 from Hartley & Zisserman
    𝐄 = 𝐊₂'*𝐅*𝐊₁
    MMatrix{3,3,Float64,3*3}(𝐄)
end
