using MultipleViewGeometry, Test, LinearAlgebra
using MultipleViewGeometry.ModuleTypes

# Construct two camera matrices and parametrise two planar surfaces.
𝐊₁ = Matrix{Float64}(I, 3, 3)
𝐑₁ = Matrix{Float64}(I, 3, 3)
𝐭₁ = [0.0, 0.0, 0.0]
𝐊₂ = Matrix{Float64}(I, 3, 3)
𝐑₂ = Matrix{Float64}(I, 3, 3)
𝐭₂ = [2.0, 2.0, 2.0]
𝐧₁ = [1.0, 0.0, 0.0]
d₁ = 10
𝐧₂ = [0.5, 0.5, 0.0]
d₂ = 15

𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# We will construct a pair of homography matrices and then construct a pair of
# projection matrices from the homographies.
𝐇₁ = construct(HomographyMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂,𝐧₁,d₁)
𝐇₂ = construct(HomographyMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂,𝐧₂,d₂)

# 1. Construct a Fundamental matrix from Camera matrices.
# 2. Construct Projection matrices from the Homography matrices
# 3. Construct a Fundamental matrix from the projection matrices.
# 4. The Fundamental matrices in step 1 and 3 should be equivalent up to sign
#    and scale.
𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

context = ProjectionMatrices(HomographyMatrices(), Chojnacki(), TwoViews())
𝐎₁, 𝐎₂ = construct(context, (𝐇₁,𝐇₂))
𝐅₂ = construct(FundamentalMatrix(), 𝐎₁, 𝐎₂)

# Ensure the matrices have the same scale and sign before comparing them.
𝐅 = 𝐅 / norm(𝐅)
𝐅 = 𝐅 / sign(𝐅[1,2])

𝐅₂ = 𝐅₂ / norm(𝐅₂)
𝐅₂ = 𝐅₂ / sign(𝐅₂[1,2])

@test 𝐅 ≈ 𝐅₂
