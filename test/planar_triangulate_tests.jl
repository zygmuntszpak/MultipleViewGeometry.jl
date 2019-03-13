using MultipleViewGeometry, Test, Random
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using LinearAlgebra
using StaticArrays

# Fix random seed.
Random.seed!(1234)
# Construct two camera matrices and parametrise two planar surfaces.
𝐊₁ = Matrix{Float64}(I, 3, 3)
𝐑₁ = Matrix{Float64}(I, 3, 3)
𝐭₁ = [-10.0, -55.0, 10.0]
𝐊₂ = Matrix{Float64}(I, 3, 3)
#𝐑₂ = Matrix{Float64}(I, 3, 3)
𝐑₂ = SMatrix{3,3,Float64,9}(rotxyz(pi/10,pi/10,pi/10))
𝐭₂ = [120.0, 120.0, 20.0]
𝐧₁ = [1.0, 0.0, 0.0]
d₁ = 25.0
𝐧₂ = [0.5, 0.5, 0.0]
d₂ = 15.0

𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# We will construct a pair of homography matrices and then construct a pair of
# projection matrices from the homographies.
𝐇₁ = construct(HomographyMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂,𝐧₁,d₁)
𝐇₂ = construct(HomographyMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂,𝐧₂,d₂)


context = ProjectionMatrices(HomographyMatrices(), Chojnacki(), TwoViews())
𝐎₁, 𝐎₂ = construct(context, (𝐇₁,𝐇₂))
𝐅 = construct(FundamentalMatrix(), 𝐎₁, 𝐎₂)

# Set of corresponding points for the first and second plane.
ℳ₁ = [Point2D(x,y) for x = 0:20:320 for y = 0:15:240]
ℳ₂ = [Point2D(x,y) for x = 320:20:640 for y = 240:15:480]

ℳ₁ʹ = similar(ℳ₁)
ℳ₂ʹ = similar(ℳ₂)
for n = 1:length(ℳ₁)
    𝐦 = hom(ℳ₁[n])
    𝐦ʹ = 𝐇₁*𝐦
    ℳ₁ʹ[n] = hom⁻¹(𝐦ʹ)
end

hom(ℳ₁[2])
hom⁻¹(𝐇₁*hom(ℳ₁[2]))

for n = 1:length(ℳ₂)
    𝐦 = hom(ℳ₂[n])
    𝐦ʹ = 𝐇₂*𝐦
    ℳ₂ʹ[n] = hom⁻¹(𝐦ʹ)
end

for n = 1:length(ℳ₁)
    m₁ = ℳ₁[n]
    m₁ʹ = ℳ₁ʹ[n]
    #Base.display(m₁,m₂)
    @show m₁,m₁ʹ
end

# 3D points corresponding to the first and second planar surface
𝒴₁ = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ₁,ℳ₁ʹ))
𝒴₂ = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ₂,ℳ₂ʹ))

N = length(𝒴₁)
for n = 1:N
    X = 𝒴₁[n]
    #@test isapprox(dot(𝐧₁,X) + d₁, 0.0; atol = 1e-12)
    Base.display(dot(𝐧₁,X) + d₁)
end

# Triangulating with the same projection matrices that were used to construct
# (ℳ,ℳʹ) should yield the same 3D points as the original 𝒳.
N = length(𝒴)
for n = 1:N
    @test  isapprox(sum(abs.(𝒳[n]-𝒴[n])/3), 0.0; atol = 1e-12)
end


#𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

# To triangulate the corresponding points using the Fundamental matrix, we first
# have to factorise the Fundamental matrix into a pair of Camera matrices. Due
# to projective ambiguity, the camera matrices are not unique, and so the
# triangulated 3D points will most probably not match the original 3D points.
# However, when working with noiseless data, the projections of the triangulated
# points should satisfy the epipolar constraint. We can use this fact to
# validate that the triangulation is correctly implemented.
𝒴 = triangulate(DirectLinearTransform(),𝐅,(ℳ₁,ℳ₁ʹ))

𝐐₁, 𝐐₂ = construct(ProjectionMatrix(),𝐅)
𝒪 = project(Pinhole(),𝐐₁,𝒴₁)
𝒪ʹ= project(Pinhole(),𝐐₂,𝒴₁)
N = length(𝒪)
for n = 1:N
    𝐦 = hom(𝒪[n])
    𝐦ʹ = hom(𝒪ʹ[n])
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-14)
end
