using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleCostFunction
using MultipleViewGeometry.ModuleTypes
using BenchmarkTools, Compat
using StaticArrays

𝒳 = [Point3DH(x,y,z,1.0)
                        for x=-1000:5:1000 for y=-1000:5:1000 for z=1:-5:-1000]

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

@time ℳ = project(Pinhole(),𝐏₁,𝒳)
@time ℳʹ = project(Pinhole(),𝐏₂,𝒳)

#@btime  project(Pinhole(),𝐏₂,𝒳)

#@time ℳ = project(Pinhole(),SMatrix{3,4}(𝐏₁),𝒳)
#@time ℳʹ = project(Pinhole(),SMatrix{3,4}(𝐏₂),𝒳)

#@btime project(Pinhole(),SMatrix{3,4}(𝐏₂),𝒳)
