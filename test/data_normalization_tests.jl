using MultipleViewGeometry, Test
using StaticArrays, LinearAlgebra

# Tests for a set of two-dimensional Cartesian points represented by homogeneous
# coordinates.
ℳ = map(Point2DH,
        [(-10.0, -10.0, 1.0),
         (-10.0,  10.0, 1.0),
         ( 10.0, -10.0, 1.0),
         ( 10.0,  10.0, 1.0)])

ℳʹ, 𝐓 = hartley_normalization(ℳ)
@test ℳʹ == map(Point2DH,
                                        [(-1.0,-1.0, 1.0),
                                         (-1.0, 1.0, 1.0),
                                         (1.0, -1.0, 1.0),
                                         (1.0,  1.0, 1.0)])
 @test 𝐓 == [0.1 0.0 -0.0;
             0.0 0.1 -0.0;
             0.0 0.0 1.0]
 @test hartley_transformation(ℳ) == [0.1 0.0 -0.0;
                                         0.0 0.1 -0.0;
                                         0.0 0.0 1.0]


# Tests for a set of three-dimensional Cartesian points represented by homogeneous
# coordinates.
ℳ = map(Point3DH,
           [(-10.0, -10.0, -10.0, 1.0),
            (-10.0, -10.0,  10.0, 1.0),
            (-10.0,  10.0, -10.0, 1.0),
            (-10.0,  10.0,  10.0, 1.0),
            ( 10.0, -10.0, -10.0, 1.0),
            ( 10.0, -10.0,  10.0, 1.0),
            ( 10.0,  10.0, -10.0, 1.0),
            ( 10.0,  10.0,  10.0, 1.0)])

ℳʹ, 𝐓 = hartley_normalization(ℳ)
@test ℳʹ == map(Point3DH,
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
@test hartley_transformation(ℳ) == [0.1 0.0 0.0 -0.0;
                                        0.0 0.1 0.0 -0.0;
                                        0.0 0.0 0.1 -0.0;
                                        0.0 0.0 0.0 1.0]
