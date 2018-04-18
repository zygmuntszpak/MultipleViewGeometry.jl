function monte_carlo_covariance(𝒞::Tuple{AbstractArray, Vararg{AbstractArray}}, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}},s,ntrial)
    ℳ, ℳʹ = collect(𝒟)
    Λ₁, Λ₂ = collect(𝒞)

    ℱ = zeros(9,ntrial)
    for itrial = 1:ntrial
      # Monte-carlo estimate of the covariance matrix.

      𝒪 = map(ℳ) do 𝐦
          Point2DH(𝑛(𝐦) + Point2DH(s*vcat(randn(2,1),0.0)))
      end
      𝒪ʹ = map(ℳʹ) do 𝐦ʹ
          Point2DH(𝑛(𝐦ʹ) + Point2DH(s*vcat(randn(2,1),0.0)))
      end

      𝐅₀ = estimate(FundamentalMatrix(),DirectLinearTransform(), (𝒪, 𝒪ʹ))
      𝐅 = estimate(FundamentalMatrix(),
                              FundamentalNumericalScheme(reshape(𝐅₀,9,1), 5, 1e-10),
                                                                (Λ₁,Λ₂), (𝒪, 𝒪ʹ))
      𝐟 = reshape(𝐅,9,1)
      𝐟 = 𝐟 / norm(𝐟)
      𝐟 = 𝐟 / sign(𝐟[2])
      ℱ[:,itrial] = 𝐟
      # for correspondence in zip(M, Mʹ)
      #     m , mʹ = correspondence
      #     # Add zero-mean Gaussian noise to the coordinates.
      #     𝐦  = 𝑛(collect(Float64,m.coords)) + s*vcat(rand(2,1),0.0)
      #     𝐦ʹ = 𝑛(collect(Float64,mʹ.coords)) + s*vcat(rand(2,1),0.0)
      #     𝒪[i] =  HomogeneousPoint(tuple(𝐦...))
      #     𝒪ʹ[i] =  HomogeneousPoint(tuple(𝐦ʹ...))
      #
      #     #@show 𝐦
      #     #@show 𝐦ʹ
      #     # Test the Fundamental Numerical Scheme on the Fundamental matrix problem.
      #     𝐅₀ = estimate(FundamentalMatrix(),DirectLinearTransform(), 𝒪, 𝒪ʹ)
      #     # 𝐅 = estimate(FundamentalMatrix(),
      #     #                         FundamentalNumericalScheme(reshape(𝐅₀,9,1), 5, 1e-10),
      #     #                         [eye(4) for i = 1:length(ℳ)],  𝒪, 𝒪ʹ)
      #     # 𝐟 = reshape(𝐅,9,1)
      #     # 𝐟 = 𝐟 / norm(𝐟)
      #     i = i + 1
      #     # 𝐟 = 𝐟 / sign(𝐟[2])
      # end
      #ℱ[:,itrial] = 𝐟
    end
    ℱ
end

#
