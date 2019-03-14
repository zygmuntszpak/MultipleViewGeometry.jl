using MultipleViewGeometry, Test, LinearAlgebra
using MultipleViewGeometry.ModuleTypes
using StaticArrays, Calculus, GeometryTypes
using MultipleViewGeometry.ModuleSyntheticData
using Random

# Fix random seed.
Random.seed!(1234)
# Construct two camera matrices and parametrise two planar surfaces.
f = 50
image_width = 640 / 10
image_height = 480 / 10
𝐊₁ = @SMatrix [f 0 0 ;
               0 f 0 ;
               0 0 1 ]
𝐑₁ = SMatrix{3,3,Float64,9}(rotxyz(0, 25*(pi/180), 0))
𝐭₁ = [-30.0, 0.0, -5.0]

𝐊₂ = @SMatrix [f 0 0 ;
               0 f 0 ;
               0 0 1 ]

𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(0, -25*(pi/180), 0))
𝐭₂ = [30.0, 0.0, 5.0]

# Normals and distance from origin
𝐧₁ = [0.0, 0.0, 1.0]
d₁ = 55.0

# Normals and distance from origin
𝐧₁ = [0.0, 0.0, 1.0]
d₁ = 55.0

𝐧₂ = [0.5, -0.2, 2.0]
d₂ = 145.0

𝒳₁ = generate_planar_points(𝐧₁,d₁, 20, 5)
𝒳₂ = generate_planar_points(𝐧₂,d₂, 20, 5)


world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
picture_basis = (Point(0.0, 0.0), Vec(-1.0, 0.0), Vec(0.0, -1.0))

camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
relocate!(camera₁, 𝐑₁, 𝐭₁)
relocate!(camera₂, 𝐑₂, 𝐭₂)

𝐑₁′, 𝐭₁′ = ascertain_pose(camera₁, world_basis... )
𝐊₁′ = obtain_intrinsics(camera₁, CartesianSystem())
𝐑₂′, 𝐭₂′ = ascertain_pose(camera₂, world_basis... )
𝐊₂′ = obtain_intrinsics(camera₂, CartesianSystem())

𝐏₁ = construct(ProjectionMatrix(),𝐊₁′,𝐑₁′,𝐭₁′)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂′,𝐑₂′,𝐭₂′)

# Set of corresponding points.
ℳ₁ = project(camera₁,𝐏₁,𝒳₁)
ℳ₁ʹ= project(camera₂,𝐏₂,𝒳₁)
ℳ₂ = project(camera₁,𝐏₁,𝒳₂)
ℳ₂ʹ= project(camera₂,𝐏₂,𝒳₂)

𝐇₁ = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ₁, ℳ₁ʹ))
𝐇₂ = estimate(HomographyMatrix(), DirectLinearTransform(), (ℳ₂, ℳ₂ʹ))


for i = zip(ℳ₁,ℳ₁ʹ)
    m, mʹ =  i
    𝐦  = hom(m)
    𝐦ʹ = hom(mʹ)
    residual = vec2antisym(𝐦ʹ)*𝐇₁*𝐦
    @test isapprox(sum(residual), 0.0; atol = 1e-7)
end

for i = zip(ℳ₂,ℳ₂ʹ)
    m, mʹ =  i
    𝐦  = hom(m)
    𝐦ʹ = hom(mʹ)
    residual = vec2antisym(𝐦ʹ)*𝐇₂*𝐦
    @test isapprox(sum(residual), 0.0; atol = 1e-7)
end


# Test the Fundamental Numerical Scheme on the Fundamental matrix problem.
Λ₁ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ₁)]
Λ₂ =  [SMatrix{3,3}(Matrix(Diagonal([1.0,1.0,0.0]))) for i = 1:length(ℳ₁ʹ)]
# 𝐇₁ʹ = estimate(HomographyMatrix(), FundamentalNumericalScheme(vec(𝐇₁), 5, 1e-10), (Λ₁,Λ₂), (ℳ₁, ℳ₁ʹ))
# 𝐇₂ʹ = estimate(HomographyMatrix(), FundamentalNumericalScheme(vec(𝐇₂), 5, 1e-10), (Λ₁,Λ₂), (ℳ₂, ℳ₂ʹ))
𝐇₁ʹ = estimate(HomographyMatrix(), FundamentalNumericalScheme(ManualEstimation(𝐇₁), 5, 1e-10), (Λ₁,Λ₂), (ℳ₁, ℳ₁ʹ))
𝐇₂ʹ = estimate(HomographyMatrix(), FundamentalNumericalScheme(ManualEstimation(𝐇₂), 5, 1e-10), (Λ₁,Λ₂), (ℳ₂, ℳ₂ʹ))

# Ensure the estimated and true matrix have the same scale and sign.
𝐇₁ = 𝐇₁ / norm(𝐇₁)
𝐇₁ = 𝐇₁ / sign(𝐇₁[1,1])

𝐇₁ʹ = 𝐇₁ʹ / norm(𝐇₁ʹ)
𝐇₁ʹ = 𝐇₁ʹ / sign(𝐇₁ʹ[1,1])

𝐇₂ = 𝐇₂ / norm(𝐇₂)
𝐇₂ = 𝐇₂ / sign(𝐇₂[1,1])

𝐇₂ʹ = 𝐇₂ʹ / norm(𝐇₂ʹ)
𝐇₂ʹ = 𝐇₂ʹ / sign(𝐇₂ʹ[1,1])

@test 𝐇₁ ≈ 𝐇₁ʹ
@test 𝐇₂ ≈ 𝐇₂ʹ

# 𝐇₁ʹʹ, fit = estimate(HomographyMatrix(), BundleAdjustment(vec(𝐇₁), 5, 1e-10), (ℳ₁, ℳ₁ʹ))
# 𝐇₂ʹʹ, fit = estimate(HomographyMatrix(), BundleAdjustment(vec(𝐇₂), 5, 1e-10), (ℳ₂, ℳ₂ʹ))
𝐇₁ʹʹ = estimate(HomographyMatrix(), BundleAdjustment(ManualEstimation(𝐇₁), 5, 1e-10), (ℳ₁, ℳ₁ʹ))
𝐇₂ʹʹ = estimate(HomographyMatrix(), BundleAdjustment(ManualEstimation(𝐇₂), 5, 1e-10), (ℳ₂, ℳ₂ʹ))

@test 𝐇₁ ≈ 𝐇₁ʹʹ
@test 𝐇₂ ≈ 𝐇₂ʹʹ
