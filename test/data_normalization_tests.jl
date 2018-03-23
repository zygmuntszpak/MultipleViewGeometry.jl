using MultipleViewGeometry, Base.Test


# Tests for a set of two-dimensional Cartesian points represented by homogeneous
# coordinates.
pts2D = map(HomogeneousPoint,
        [(-10.0, -10.0, 1.0),
         (-10.0,  10.0, 1.0),
         ( 10.0, -10.0, 1.0),
         ( 10.0,  10.0, 1.0)])

pts2Dʹ, 𝐓 = hartley_normalization(pts2D)
@test pts2Dʹ == map(HomogeneousPoint,
                                        [(-1.0,-1.0, 1.0),
                                         (-1.0, 1.0, 1.0),
                                         (1.0, -1.0, 1.0),
                                         (1.0,  1.0, 1.0)])
 @test 𝐓 == [0.1 0.0 -0.0;
             0.0 0.1 -0.0;
             0.0 0.0 1.0]
 @test hartley_transformation(pts2D) == [0.1 0.0 -0.0;
                                         0.0 0.1 -0.0;
                                         0.0 0.0 1.0]


# Tests for a set of three-dimensional Cartesian points represented by homogeneous
# coordinates.
pts3D = map(HomogeneousPoint,
           [(-10.0, -10.0, -10.0, 1.0),
            (-10.0, -10.0,  10.0, 1.0),
            (-10.0,  10.0, -10.0, 1.0),
            (-10.0,  10.0,  10.0, 1.0),
            ( 10.0, -10.0, -10.0, 1.0),
            ( 10.0, -10.0,  10.0, 1.0),
            ( 10.0,  10.0, -10.0, 1.0),
            ( 10.0,  10.0,  10.0, 1.0)])

pts3Dʹ, 𝐓 = hartley_normalization(pts3D)
@test pts3Dʹ == map(HomogeneousPoint,
                                         [(-1.0,-1.0, -1.0, 1.0),
                                          (-1.0,-1.0,  1.0, 1.0),
                                          (-1.0, 1.0, -1.0, 1.0),
                                          (-1.0, 1.0,  1.0, 1.0),
                                          (1.0, -1.0, -1.0, 1.0),
                                          (1.0, -1.0,  1.0, 1.0),
                                          (1.0,  1.0, -1.0, 1.0),
                                          (1.0,  1.0,  1.0, 1.0)])
@test 𝐓 == [0.1 0.0 0.0 -0.0;
            0.0 0.1 0.0 -0.0;
            0.0 0.0 0.1 -0.0;
            0.0 0.0 0.0 1.0]
@test hartley_transformation(pts3D) == [0.1 0.0 0.0 -0.0;
                                        0.0 0.1 0.0 -0.0;
                                        0.0 0.0 0.1 -0.0;
                                        0.0 0.0 0.0 1.0]
