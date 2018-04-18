function moments(entity::FundamentalMatrix, 𝒟::Tuple{AbstractArray, Vararg{AbstractArray}})
    # ʹ : CTRL + SHIFT + 02b9
    ℳ, ℳʹ = collect(𝒟)
    N = length(ℳ)
    if (N != length(ℳʹ))
           throw(ArgumentError("There should be an equal number of points for each view."))
    end
    𝐀 =  @SMatrix zeros(9,9)
    for n = 1:N
        𝐦  = 𝑛(ℳ[n])
        𝐦ʹ = 𝑛(ℳʹ[n])
        𝐀 = 𝐀 + (𝐦*𝐦') ⊗ (𝐦ʹ*𝐦ʹ')
    end
    𝐀/N
end
