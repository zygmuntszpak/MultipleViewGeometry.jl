using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using BenchmarkTools, Compat
using StaticArrays

# Fix random seed.
srand(1234)

𝒳 = [Point3DH(x,y,z,1.0)
                        for x=-1:0.5:10 for y=-1:0.5:10 for z=2:-0.1:1]

# Intrinsic and extrinsic parameters of camera one.
𝐊₁ = @SMatrix eye(3)
𝐑₁ = @SMatrix eye(3)
𝐭₁ =  @SVector [0.0, 0.0, -10]

# Intrinsic and extrinsic parameters of camera two.
𝐊₂ = @SMatrix eye(3)
𝐑₂ = @SMatrix eye(3) #SMatrix{3,3,Float64,9}(rotxyz(pi/10,pi/10,pi/10))
𝐭₂ = @SVector [10.0, 10.0, -10.0]

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁,𝐑₁,𝐭₁)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂,𝐑₂,𝐭₂)

# Set of corresponding points.
ℳ = project(Pinhole(),𝐏₁,𝒳)
ℳʹ = project(Pinhole(),𝐏₂,𝒳)

𝒴 = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ,ℳʹ))

# Triangulating with the same projection matrices that were used to construct
# (ℳ,ℳʹ) should yield the same 3D points as the original 𝒳.
N = length(𝒴)
for n = 1:N
    @test  isapprox(sum(abs.(𝒳[n]-𝒴[n])/4), 0.0; atol = 1e-12)
end


𝐅 = construct(FundamentalMatrix(),𝐊₁,𝐑₁,𝐭₁,𝐊₂,𝐑₂,𝐭₂)

# To triangulate the corresponding points using the Fundamental matrix, we first
# have to factorise the Fundamental matrix into a pair of Camera matrices. Due
# to projective ambiguity, the camera matrices are not unique, and so the
# triangulated 3D points will most probably not match the original 3D points.
# However, when working with noiseless data, the projections of the triangulated
# points should satisfy the epipolar constraint. We can use this fact to
# validate that the triangulation is correctly implemented.
𝒴 = triangulate(DirectLinearTransform(),𝐅,(ℳ,ℳʹ))

𝐐₁, 𝐐₂ = construct(ProjectionMatrix(),𝐅)
𝒪 = project(Pinhole(),𝐐₁,𝒴)
𝒪ʹ= project(Pinhole(),𝐐₂,𝒴)
N = length(𝒪)
for n = 1:N
    𝐦 = 𝒪[n]
    𝐦ʹ = 𝒪ʹ[n]
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-14)
end
