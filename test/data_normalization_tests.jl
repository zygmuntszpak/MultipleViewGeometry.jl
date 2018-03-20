using MultipleViewGeometry, Base.Test


# Tests for a set of two-dimensional Cartesian points represented by homogeneous
# coordinates.
pts2D = map(HomogeneousPoint,
        [(-10.0, -10.0, 1.0),
         (-10.0,  10.0, 1.0),
         ( 10.0, -10.0, 1.0),
         ( 10.0,  10.0, 1.0)])

@test hartley_normalization(pts2D) == map(HomogeneousPoint,
                                        [(-1.0,-1.0, 1.0),
                                         (-1.0, 1.0, 1.0),
                                         (1.0, -1.0, 1.0),
                                         (1.0,  1.0, 1.0)])

 @test hartley_transformation(pts2D) == eye(3)


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

@test hartley_normalization(pts3D) == map(HomogeneousPoint,
                                         [(-1.0,-1.0, -1.0, 1.0),
                                          (-1.0,-1.0,  1.0, 1.0),
                                          (-1.0, 1.0, -1.0, 1.0),
                                          (-1.0, 1.0,  1.0, 1.0),
                                          (1.0, -1.0, -1.0, 1.0),
                                          (1.0, -1.0,  1.0, 1.0),
                                          (1.0,  1.0, -1.0, 1.0),
                                          (1.0,  1.0,  1.0, 1.0)])

@test hartley_transformation(pts3D) == eye(4)
