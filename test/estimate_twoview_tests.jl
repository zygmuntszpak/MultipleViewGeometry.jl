using MultipleViewGeometry, Base.Test
using StaticArrays
# Tests for fundamental matrix estimation


𝒳 = [Point3DH(x,y,z,1.0)
                        for x=-100:5:100 for y=-100:5:100 for z=1:-50:-100]

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
𝐅 = estimate(FundamentalMatrix(), DirectLinearTransform(), (ℳ, ℳʹ))
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
    𝐦  = 𝑛(m)
    𝐦ʹ = 𝑛(mʹ)
    residual[i] = (𝐦ʹ'*𝐅*𝐦) 
end

@test isapprox(sum(residual), 0.0; atol = 1e-7)

# Test the Fundamental Numerical Scheme on the Fundamental matrix problem.
Λ₁ =  [SMatrix{3,3}(diagm([1.0,1.0,0.0])) for i = 1:length(ℳ)]sum(residual)
Λ₂ =  [SMatrix{3,3}(diagm([1.0,1.0,0.0])) for i = 1:length(ℳ)]
𝐅₀ = estimate(FundamentalMatrix(),DirectLinearTransform(),  (ℳ, ℳʹ))
𝐅 = estimate(FundamentalMatrix(),
                        FundamentalNumericalScheme(reshape(𝐅₀,9,1), 5, 1e-10),
                                                          (Λ₁,Λ₂), (ℳ, ℳʹ))

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
