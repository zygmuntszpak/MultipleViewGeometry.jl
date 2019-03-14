function estimate(entity::HomographyMatrix, method::DirectLinearTransform, 𝓓::Tuple{Vector{Vector{T₁}} where T₁ <: AbstractArray, Vector{Vector{T₂}} where T₂ <: AbstractArray})
    𝓜, 𝓜ʹ =  𝓓
    𝓡 = Vector{SMatrix{3,3,Float64,9}}(undef,length(𝓓))
    for k = 1:length(𝓓)
        𝓡[k] = estimate(entity, method, (𝓜[k], 𝓜ʹ[k]))
    end
    𝓡
end

function estimate(entity::HomographyMatrix, method::DirectLinearTransform, 𝒟::Tuple{AbstractArray, AbstractArray})
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

function estimate(entity::HomographyMatrix, method::FundamentalNumericalScheme,  𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, AbstractArray})
    ℳ, ℳʹ = 𝒟
    Λ₀, Λ₀ʹ = 𝒞
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end
    if (N != length(Λ₀) || N != length(Λ₀ʹ) )
          throw(ArgumentError("There should be a covariance matrix for each point correspondence."))
    end

    # Initial estimate which will be used to seed the fundmamental numerical scheme.
    𝛉₀ = vec(estimate(HomographyMatrix(), method.seed, 𝒟))

    # Map corresponding points to the normalized coordinate system.
    𝒪, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
    𝒪ʹ, 𝐓ʹ = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳʹ)
    # Map seed to the normalized coordinate system.
    𝛉 = (inv(𝐓') ⊗ 𝐓ʹ) * 𝛉₀
    𝛉 = 𝛉 / norm(𝛉)

    # Map covariance matrices to the normalized coordinate system.
    Λ₁ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ , 𝐓)
    Λ₁ʹ = transform(CovarianceMatrices(), CanonicalToHartley(), Λ₀ʹ , 𝐓ʹ)
    for i = 1:method.max_iter
        𝐗 = X(AML(),HomographyMatrix(), 𝛉, (Λ₁,Λ₁ʹ), (𝒪, 𝒪ʹ))
        λ, 𝛉⁺ = smallest_eigenpair(Symmetric(𝐗/N))
        𝛉 = reshape(𝛉⁺, length(𝛉⁺), 1)
    end
    𝐡 = (𝐓' ⊗ inv(𝐓ʹ)) * 𝛉
    𝐇 = reshape(𝐡,(3,3))
end



function estimate(entity::HomographyMatrix, method::BundleAdjustment,  𝒟::Tuple{AbstractArray, AbstractArray})
    ℳ, ℳʹ = 𝒟
    N = length(ℳ)
    if (N != length(ℳʹ))
          throw(ArgumentError("There should be an equal number of points for each view."))
    end

    # Initial estimate which will be used to seed the fundmamental numerical scheme.
    𝐇 = estimate(HomographyMatrix(), method.seed, 𝒟)

    #𝐇 = SMatrix{3,3,Float64,9}(reshape(method.𝛉₀,(3,3)))
    𝐈 = SMatrix{3,3}(1.0I)
    # Construct a length-(9+2*N) vector consisting of the homography matrix
    # (the first nine dimensions), as well as N two-dimensional points in the
    # first view (the remaining dimensions).
    𝛉 = pack(HomographyMatrix(), 𝐇, ℳ)

    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    pts = Matrix{Float64}(undef,4,N)
    for n = 1:N
        pts[index₁,n] = ℳ[n][index₁]
        pts[index₂,n] = ℳʹ[n][index₁]
    end
    fit = curve_fit(model_homography, jacobian_model_homography,  𝐈, reshape(reinterpret(Float64,vec(pts)),(4*N,)) , 𝛉; show_trace = false)
    #fit = curve_fit(model_homography, 𝐈, reshape(reinterpret(Float64,vec(pts)),(4*N,)) , 𝛉; show_trace = false)
    𝐇₊  = SMatrix{3,3,Float64,9}(reshape(fit.param[1:9],(3,3)))
    𝐇₊
end

function model_homography(𝐈,𝛉)
    # Nine parameters for the projection matrix, and 2 parameters per 2D point.
    N = Int((length(𝛉) - 9) / 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    reprojections = Matrix{Float64}(undef,4,N)
    𝛉v = @view 𝛉[1:9]
    𝐇 = SMatrix{3,3,Float64,9}(reshape(𝛉v,(3,3)))
    i = 10
    for n = 1:N
        # Extract 2D point and convert to homogeneous coordinates
        𝐦 = hom(SVector{2,Float64}(𝛉[i],𝛉[i+1]))
        reprojections[index₁,n] = hom⁻¹(𝐦)
        reprojections[index₂,n] = hom⁻¹(𝐇 * 𝐦)
        i = i + 2
    end
    reshape(reinterpret(Float64,vec(reprojections)),(4*N,))
end

function jacobian_model_homography(𝐈,𝛉)
    # Twelve parameters for the projection matrix, and 2 parameters per 2D point.
    N = Int((length(𝛉) - 9) / 2)
    index₁ = SVector(1,2)
    index₂ = SVector(3,4)
    reprojections = Matrix{Float64}(undef,4,N)
    𝛉v = @view 𝛉[1:9]
    𝐇 = SMatrix{3,3,Float64,9}(reshape(𝛉v,(3,3)))
    𝐉 = zeros(4*N,9+2*N)
    # Create a view of the jacobian matrix 𝐉 and reshape it so that
    # it will be more convenient to index into the appropriate entries
    # whilst looping over all of the data points.
    𝐉v = reshape(reinterpret(Float64,𝐉), 4, N, 9+2*N)
    𝐀 = SMatrix{2,3,Float64,6}(1,0,0,1,0,0)
    𝐈₃ = SMatrix{3,3}(1.0I)
    i = 10
    for n = 1:N
        # Extract 3D point and convert to homogeneous coordinates.
        𝐦 = hom(SVector{2,Float64}(𝛉[i], 𝛉[i+1]))

        # Derivative of residual in first and second image w.r.t 2D point in the
        # first image.
        ∂𝐫₁_d𝐦 = 𝐀 * 𝐈₃
        ∂𝐫₂_d𝐦 = 𝐀 * ∂hom⁻¹(𝐇 * 𝐦) * 𝐇

        # Derivative of residual in second image w.r.t homography martix.
        # ∂𝐫₁_d𝐇 is the zero vector.
        ∂𝐫₂_d𝐇 = 𝐀 * ∂hom⁻¹(𝐇  * 𝐦) * (𝐦' ⊗ 𝐈₃)

        𝐉v[index₂,n,1:9] = ∂𝐫₂_d𝐇
        𝐉v[index₁,n,i:i+1] = ∂𝐫₁_d𝐦[:,1:2]
        𝐉v[index₂,n,i:i+1] = ∂𝐫₂_d𝐦[:,1:2]
        i = i + 2
    end
    𝐉
end

# Construct a parameter vector consisting of a homography matrix and 2D points.
function pack(entity::HomographyMatrix, 𝐇::AbstractArray, ℳ::AbstractArray)
    N = length(ℳ)
    𝛉 = Vector{Float64}(undef,9 + N*2)
    𝛉[1:9] = Array(vec(𝐇))
    i = 10
    for n = 1:N
        𝛉[i:i+1] = ℳ[n][1:2]
        i = i + 2
    end
    𝛉
end
