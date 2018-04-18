using MultipleViewGeometry, Base.Test
using BenchmarkTools
using StaticArrays


ℳ = [Point2DH(x,y,1)  for x=-1000:0.5:1000 for y=-1000:0.5:1000]
#ℳ = [Point2DH(x,y,1)  for x=-10:0.5:10 for y=-10:0.5:10]

transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)
ℳʹ, 𝐓 = transform(HomogeneousCoordinates(),CanonicalToHartley(),ℳ)


# @time Λ =  [MMatrix{3,3}(diagm([1.0,1.0,0.0])) for i = 1:length(ℳ)]
# @time transform(CovarianceMatrices(), CanonicalToHartley(), Λ , tuple(𝐓))
# @time Z = transform(CovarianceMatrices(), CanonicalToHartley(), Λ , tuple(𝐓))

@time Λ =  [SMatrix{3,3}(diagm([1.0,1.0,0.0])) for i = 1:length(ℳ)]
@time transform(CovarianceMatrices(), CanonicalToHartley(), Λ , 𝐓)
@time Z = transform(CovarianceMatrices(), CanonicalToHartley(), Λ , 𝐓)

@show Λ[1]
@show Z[1]
