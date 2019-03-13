using MultipleViewGeometry, Test
using MultipleViewGeometry.ModuleTypes


ℳ = [Point2D(rand(1,2)) for i = 1:10000]

moments(FundamentalMatrix(), (ℳ,ℳ))

@benchmark moments(HomographyMatrix(), (ℳ,ℳ))

# TODO

FundamentalMatrix()
HomographyMatrix()
