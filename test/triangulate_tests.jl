using MultipleViewGeometry, Test, Random
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using MultipleViewGeometry.ModuleConstraints
using MultipleViewGeometry.ModuleConstruct
using LinearAlgebra
using StaticArrays

# Fix random seed.
Random.seed!(1234)


𝒳 = [Point3D(x,y,rand(50:100)) for x = -100:5:100 for y = -100:5:100]
𝒳 = 𝒳[1:50:end]


# Specify the coordinate systems of the world, the camera frame and the picture
# plane.
world_basis = (Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0))
camera_basis = (Point(0.0, 0.0, 0.0), Vec(-1.0, 0.0, 0.0), Vec(0.0, -1.0, 0.0), Vec(0.0, 0.0, 1.0))
picture_basis = (Point(0.0, 0.0), Vec(-1.0, 0.0), Vec(0.0, -1.0))

# The focal length for both cameras is one.
f = 1
image_width = 640 / 10
image_height = 480 / 10

camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)

# Rotate and translate camera one.
𝐑₁ = Matrix{Float64}(I,3,3)
𝐭₁ = [-50.0, -2.0, 0.0]
relocate!(camera₁, 𝐑₁, 𝐭₁)

# Rotate and translate camera two.
𝐑₂ = Matrix{Float64}(I,3,3)
𝐭₂ = [50.0, 2.0, 0.0]
relocate!(camera₂, 𝐑₂, 𝐭₂)


𝐑₁′, 𝐭₁′ = ascertain_pose(camera₁, world_basis... )
𝐊₁′ = obtain_intrinsics(camera₁, CartesianSystem())
𝐑₂′, 𝐭₂′ = ascertain_pose(camera₂, world_basis... )
𝐊₂′ = obtain_intrinsics(camera₂, CartesianSystem())

# Camera projection matrices.
𝐏₁ = construct(ProjectionMatrix(),𝐊₁′,𝐑₁′,𝐭₁′)
𝐏₂ = construct(ProjectionMatrix(),𝐊₂′,𝐑₂′,𝐭₂′)


# Set of corresponding points.
ℳ = project(camera₁,𝐏₁,𝒳)
ℳʹ = project(camera₂,𝐏₂,𝒳)

𝒴 = triangulate(DirectLinearTransform(),𝐏₁,𝐏₂,(ℳ,ℳʹ))

# Triangulating with the same projection matrices that were used to construct
# (ℳ,ℳʹ) should yield the same 3D points as the original 𝒳.
N = length(𝒴)
for n = 1:N
    @test  isapprox(sum(abs.(𝒳[n]-𝒴[n])/3), 0.0; atol = 1e-11)
end


𝐅 = construct(FundamentalMatrix(),𝐊₁′,𝐑₁′,-𝐭₁′,𝐊₂′,𝐑₂′,-𝐭₂′)

# To triangulate the corresponding points using the Fundamental matrix, we first
# have to factorise the Fundamental matrix into a pair of Camera matrices. Due
# to projective ambiguity, the camera matrices are not unique, and so the
# triangulated 3D points will most probably not match the original 3D points.
# However, when working with noiseless data, the projections of the triangulated
# points should satisfy the epipolar constraint. We can use this fact to
# validate that the triangulation is correctly implemented.
𝒴 = triangulate(DirectLinearTransform(),𝐅,(ℳ,ℳʹ))

𝐐₁, 𝐐₂ = construct(ProjectionMatrix(),𝐅)
# Because we constructed the projection matrices from the fundamental matrix we
# don't know the intrinsics or extrinsics of the camera.
# The current API requires us to construct a CameraModel type for the `project`
# function. TODO: Need to revisit this.
camera₁ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
camera₂ = Pinhole(image_width, image_height, f, camera_basis..., picture_basis...)
𝒪 = project(camera₁,𝐐₁,𝒴)
𝒪ʹ= project(camera₂,𝐐₂,𝒴)
N = length(𝒪)
for n = 1:N
    𝐦 = hom(𝒪[n])
    𝐦ʹ = hom(𝒪ʹ[n])
    @test  isapprox(𝐦'*𝐅*𝐦ʹ, 0.0; atol = 1e-12)
end
