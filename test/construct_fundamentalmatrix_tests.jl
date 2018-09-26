using MultipleViewGeometry, Test, LinearAlgebra

𝐊₁ = Matrix{Float64}(I, 3, 3)
𝐑₁ = Matrix{Float64}(I, 3, 3)
𝐭₁ = [1.0, 1.0, 1.0]
𝐊₂ = Matrix{Float64}(I, 3, 3)
𝐑₂ = Matrix{Float64}(I, 3, 3)
𝐭₂ = [2.0, 2.0, 2.0]

@test construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂) == vec2antisym([-1,-1,-1])
