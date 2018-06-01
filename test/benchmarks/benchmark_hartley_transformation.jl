using MultipleViewGeometry, Base.Test
using BenchmarkTools
using StaticArrays


ℳ = [Point2DH(x,y,1)  for x=-1000:0.5:1000 for y=-1000:0.5:1000]
𝐓 = hartley_transformation(ℳ)
@time hartley_transformation(ℳ)
@time hartley_transformation(ℳ)

@time hartley_normalization(ℳ)
@time hartley_normalization(ℳ)

𝒪, 𝐓 = hartley_normalization(ℳ)

ℳ[1]
𝒪[1]

@time transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
