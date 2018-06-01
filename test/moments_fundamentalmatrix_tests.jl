using MultipleViewGeometry, Base.Test
using MultipleViewGeometry.ModuleTypes


ℳ = [Point2DH(rand(1,3)) for i = 1:10000]

moments(FundamentalMatrix(), (ℳ,ℳ))

# TODO
