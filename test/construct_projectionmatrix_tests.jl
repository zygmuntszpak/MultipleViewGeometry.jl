using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleTypes

𝐊 = eye(3)
𝐑 = eye(3)
𝐭 = [1.0, 1.0, 1.0]

@test construct(ProjectionMatrix(),𝐊,𝐑,𝐭) == [eye(3) -ones(3)]

# 1. Construct a Fundamental matrix from Camera matrices.
# 2. Construct projection matrices from the Fundamental matrix.
# 3. Construct a Fundamental matrix from the projection matrices.
# 4. The Fundamental matrices in step 2 and 3 should be equivalent up to sign
#    and scale.
𝐊₁ = eye(3)
𝐑₁ = eye(3)
𝐭₁ = [1.0, 1.0, 1.0]
𝐊₂ = eye(3)
𝐑₂ = eye(3)
𝐭₂ = [2.0, 2.0, 2.0]
𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)
𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

𝐎₁, 𝐎₂ = construct(ProjectionMatrix(),𝐅)
𝐅₂ = construct(FundamentalMatrix(),𝐎₁, 𝐎₂)

# Ensure the matrices have the same scale and sign before comparing them.
𝐅 = 𝐅 / norm(𝐅)
𝐅 = 𝐅 / sign(𝐅[1,2])

𝐅₂ = 𝐅₂ / norm(𝐅₂)
𝐅₂ = 𝐅₂ / sign(𝐅₂[1,2])

@test 𝐅 ≈ 𝐅₂
