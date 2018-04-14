using MultipleViewGeometry, Base.Test

ℳ = map(HomogeneousPoint,
        [(-10.0, -10.0, 1.0),
         (-10.0,  10.0, 1.0),
         ( 10.0, -10.0, 1.0),
         ( 10.0,  10.0, 1.0)])

ℳ = map(HomogeneousPoint,
          [(-20.0, -20.0, 1.0),
           (-20.0,  20.0, 1.0),
           ( 20.0, -20.0, 1.0),
           ( 20.0,  20.0, 1.0)])

moments(FundamentalMatrix(), (ℳ,ℳ)...)

# TODO
