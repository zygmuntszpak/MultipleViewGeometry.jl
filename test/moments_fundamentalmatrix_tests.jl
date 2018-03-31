using MultipleViewGeometry, Base.Test

pts1 = map(HomogeneousPoint,
        [(-10.0, -10.0, 1.0),
         (-10.0,  10.0, 1.0),
         ( 10.0, -10.0, 1.0),
         ( 10.0,  10.0, 1.0)])

pts2 = map(HomogeneousPoint,
          [(-20.0, -20.0, 1.0),
           (-20.0,  20.0, 1.0),
           ( 20.0, -20.0, 1.0),
           ( 20.0,  20.0, 1.0)])

moments(FundamentalMatrix(), (pts1,pts2)...)

# TODO
