using MultipleViewGeometry, Test, LinearAlgebra

# In the case of calibrated cameras one can assume that
# 𝐊₁ = 𝐊₂ = 𝐈, 𝐭₁ = [0,0,0],  𝐑₁ = 𝐈 ,  𝐑₂ = 𝐑 and 𝐭 = -𝐑𝐭₂ so that the
# direct nRt representation 𝐇 = -d𝐑 + 𝐭𝐧' holds.
# We will verify that our general function for constructing a homography
# produces the same result as the direct nRt representation for calibrated
# cameras.
𝐊₁ = Matrix{Float64}(I, 3, 3)
𝐑₁ = Matrix{Float64}(I, 3, 3)
𝐭₁ = [0.0, 0.0, 0.0]
𝐊₂ = Matrix{Float64}(I, 3, 3)
𝐑₂ = Matrix{Float64}(I, 3, 3)
𝐭₂ = [2.0, 2.0, 2.0]
𝐧 = [1.0, 0.0, 0.0]
d = 10

𝐑 = 𝐑₂
𝐭 = -𝐑*𝐭₂
𝐇₀ = -d*𝐑 + 𝐭*𝐧'

𝐇 = construct(HomographyMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂,𝐧,d)

# Homographies are equivalent up to scale and sign so we need to normalise them
# before comparing them.
𝐇₀ = 𝐇₀ / norm(𝐇₀)
𝐇₀ = 𝐇₀ / sign(𝐇₀[3,3])

𝐇 = 𝐇 / norm(𝐇)
𝐇 = 𝐇 / sign(𝐇[3,3])


@test all(𝐇 .≈ 𝐇₀)
