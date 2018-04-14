using MultipleViewGeometry, Base.Test

# Tests for fundamental matrix estimation

# A rectangular array of 3D points represented in homogeneous coordinates
𝒳 = [HomogeneousPoint(Float64.((x,y,z,1.0),RoundDown))
                        for x=-100:10:100 for y=-100:10:100 for z=1:-100:-1000]

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

# Estimate of the fundamental matrix and the true fundamental matrix.
𝐅 = estimate(FundamentalMatrix(), DirectLinearTransform(), ℳ, ℳʹ)
𝐅ₜ = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

# Ensure the estimated and true matrix have the same scale and sign.
𝐅 = 𝐅 / norm(𝐅)
𝐅 = 𝐅 / sign(𝐅[1,2])
𝐅ₜ = 𝐅ₜ / norm(𝐅ₜ)
𝐅ₜ = 𝐅ₜ / sign(𝐅ₜ[1,2])

@test 𝐅 ≈ 𝐅ₜ

# Check that the fundamental matrix satisfies the corresponding point equation.
npts = length(ℳ)
residual = zeros(Float64,npts,1)
for correspondence in zip(1:length(ℳ),ℳ, ℳʹ)
    i, m , mʹ = correspondence
    𝐦  = 𝑛(collect(Float64,m.coords))
    𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))
    residual[i] = 𝐦ʹ'*𝐅*𝐦
end

@test isapprox(sum(residual), 0.0; atol = 1e-9)

# Test the Fundamental Numerical Scheme on the Fundamental matrix problem.
𝐅₀ = estimate(FundamentalMatrix(),DirectLinearTransform(), ℳ, ℳʹ)
𝐅 = estimate(FundamentalMatrix(),
                        FundamentalNumericalScheme(reshape(𝐅₀,9,1), 5, 1e-10),
                        [eye(4) for i = 1:length(ℳ)], ℳ, ℳʹ)

𝐅ₜ = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)
# Ensure the estimated and true matrix have the same scale and sign.
𝐅 = 𝐅 / norm(𝐅)
𝐅 = 𝐅 / sign(𝐅[1,2])
𝐅ₜ = 𝐅ₜ / norm(𝐅ₜ)
𝐅ₜ = 𝐅ₜ / sign(𝐅ₜ[1,2])

@test 𝐅 ≈ 𝐅ₜ

# The way the Taubin estimate is implemented is numerically unstable
# for noiseless data.

# # Estimate of the fundamental matrix and the true fundamental matrix.
# 𝐅 = estimate(FundamentalMatrix(),Taubin(), ℳ, ℳʹ)
# 𝐅ₜ = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)
#
# # Ensure the estimated and true matrix have the same scale and sign.
# 𝐅 = 𝐅 / norm(𝐅)
# 𝐅 = 𝐅 / sign(𝐅[1,2])
# 𝐅ₜ = 𝐅ₜ / norm(𝐅ₜ)
# 𝐅ₜ = 𝐅ₜ / sign(𝐅ₜ[1,2])
#
# @test 𝐅 ≈ 𝐅ₜ


# m = ℳ[1]
# mʹ = ℳʹ[1]
# 𝐦  = 𝑛(collect(Float64,m.coords))
# 𝐦ʹ = 𝑛(collect(Float64,mʹ.coords))

# 𝐞₁ = [1.0 0.0 0.0]'
# 𝐞₂ = [0.0 1.0 0.0]'
# (𝐦*𝐦') ⊗ (𝐦ʹ*𝐦ʹ')
# z = [kron(𝐞₁, 𝐦ʹ)]
# z = [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
# ∂ₓ𝐮 = [(𝐞₁ ⊗ 𝐦ʹ) (𝐞₂ ⊗ 𝐦ʹ) (𝐦 ⊗ 𝐞₁) (𝐦 ⊗ 𝐞₂)]
