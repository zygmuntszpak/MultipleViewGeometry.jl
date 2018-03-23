using MultipleViewGeometry, Base.Test

# Tests for fundamental matrix estimation
ℳ = map(HomogeneousPoint, [(-10.0, -10.0, 1.0),
                               (-10.0,  10.0, 1.0),
                               ( 10.0,-10.0, 1.0),
                               ( 10.0,  10.0, 1.0)])

ℳʹ = map(HomogeneousPoint, [(10.0, -10.0, 1.0),
                               (10.0,  10.0, 1.0),
                               ( 20.0,-10.0, 1.0),
                               ( 20.0,  10.0, 1.0)])

F = estimate(FundamentalMatrix(), ℳ, ℳʹ)
