using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleCostFunction

# Test for cost functions.

# Test cost function on Fundamental matrix estimation.

# A rectangular array of 3D points represented in homogeneous coordinates
# 𝒳 = [HomogeneousPoint(Float64.((x,y,z,1.0),RoundDown))
#                         for x=-100:10:100 for y=-100:10:100 for z=1:-100:-1000]

𝒳 = [HomogeneousPoint(Float64.((x,y,z,1.0),RoundDown))
                        for x=-100:50:100 for y=-100:50:100 for z=1:-500:-1000]

# Intrinsic and extrinsic parameters of camera one.
𝐊₁ = eye(3)
𝐑₁ = eye(3)
𝐭₁ = [0.0, 0.0, 0.0]

# Intrinsic and extrinsic parameters of camera two.
𝐊₂ = eye(3)
𝐑₂ = eye(3)
𝐭₂ = [100.0, 2.0, -100.0]

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# Set of corresponding points.
ℳ = project(Pinhole(),𝐏₁,𝒳)
ℳʹ = project(Pinhole(),𝐏₂,𝒳)


𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

# Ensure the estimated and true matrix have the same scale and sign.
𝐅 = 𝐅 / norm(𝐅)
𝐅 = 𝐅 / sign(𝐅[1,2])
𝐟 = reshape(𝐅,9,1)
Jₐₘₗ =  cost(AML(),FundamentalMatrix(), 𝐟 ,
                                         [eye(4) for i = 1:length(ℳ)], ℳ, ℳʹ)

@test isapprox(Jₐₘₗ, 0.0; atol = 1e-14)

# Verify that the the vectorised fundamental matrix is in the null space of X
𝐗 = X(AML(),FundamentalMatrix(), reshape(𝐅,9,1),
                                         [eye(4) for i = 1:length(ℳ)], ℳ, ℳʹ)

# The true parameters should lie in the null space of the matrix X.
@test all(isapprox.(𝐗 * 𝐟, 0.0; atol = 1e-10))

# matches = ℳ, ℳʹ
# Λ = [eye(4) for i = 1:length(ℳ)]
# zip(matches, Λ)
#
# 𝐗*𝐟
