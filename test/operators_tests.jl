using MultipleViewGeometry, Base.Test
using StaticArrays

# Vectors are scaled so that the last component is unity.
v = [2,2,2]
@test 𝑛(v) == [1.0, 1.0, 1.0]

v = Point2DH([2,2,2])
@test 𝑛(v) == Point2DH([1.0, 1.0, 1.0])

# Vectors which represent points at infinity are unchanged.
v = [2,2,0]
@test 𝑛(v) == [2,2,0]

v = Point2DH([2,2,0])
@test 𝑛(v) == Point2DH([2,2,0])

v = [1, 2, 3]
@test vec2antisym(v) == @SMatrix [ 0 -3  2;
                                   3  0 -1;
                                  -2  1  0]

 v = Point2DH([1, 2, 3])
 @test vec2antisym(v) == @SMatrix [ 0 -3  2;
                                    3  0 -1;
                                   -2  1  0]
