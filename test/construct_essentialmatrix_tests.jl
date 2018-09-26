using MultipleViewGeometry, MultipleViewGeometry.ModuleRotation, LinearAlgebra, Test

# Intrinsic and extrinsic parameters for the first camera.
𝐊₁ = zeros(3,3)
𝐊₁[1,1] = 10
𝐊₁[2,2] = 10
𝐊₁[3,3] = 1
𝐑₁ = rotxyz(deg2rad(10), deg2rad(15), deg2rad(45))
𝐭₁ = [-250.0, 0.0, 2500.0]

# Intrinsic and extrinsic parameters for the second camera.
𝐊₂ = zeros(3,3)
𝐊₂[1,1] = 5
𝐊₂[2,2] = 5
𝐊₂[3,3] = 1
𝐑₂ =  rotxyz(deg2rad(10), deg2rad(15), deg2rad(45))
𝐭₂ = [250.0,   0.0, 2500.0]

𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)
𝐄 = construct(EssentialMatrix(),𝐅, 𝐊₁, 𝐊₂)

# Result 9.17 of R. Hartley and A. Zisserman, “Two-View Geometry,” Multiple View Geometry in Computer Vision
# A 3 by 3 matrix is an essential matrix if and only if two of its singular values
# are equal, and the third is zero.
U, S , V = svd(𝐄)

@test isapprox.(S[1], S[2]; atol = 1e-14)
@test isapprox.(S[3], 0.0; atol = 1e-10)
