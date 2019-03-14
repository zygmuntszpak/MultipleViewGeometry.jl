function estimate(entity::HomographyMatrix, method::DirectLinearTransform, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    ℳ, ℳʹ =  𝒟
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    𝒪, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    𝒪ʹ, 𝐓ʹ  = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    𝐀 = moments(HomographyMatrix(), (𝒪, 𝒪ʹ))
    λ, h = smallest_eigenpair(Symmetric(𝐀))
    𝐇 = reshape(h,(3,3))
    𝐇 = SMatrix{3,3,Float64,9}(𝐇 / norm(𝐇))
    # Transform estimate back to the original (unnormalised) coordinate system.
    inv(𝐓ʹ)*𝐇*𝐓
end
