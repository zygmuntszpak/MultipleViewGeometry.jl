using MultipleViewGeometry, Test, LinearAlgebra

@test rotx(0.0) == Matrix(1.0I, 3, 3)
@test rotx(2.0*pi) ≈  Matrix(1.0I, 3, 3)

@test roty(0.0) == Matrix(1.0I, 3, 3)
@test roty(2.0*pi) ≈  Matrix(1.0I, 3, 3)

@test rotz(0.0) == Matrix(1.0I, 3, 3)
@test rotz(2.0*pi) ≈  Matrix(1.0I, 3, 3)

@test rotxyz(0.0, 0.0, 0.0) == Matrix(1.0I, 3, 3)
@test rotxyz(2.0*pi, 2.0*pi, 2.0*pi) ≈  Matrix(1.0I, 3, 3)

@test rodrigues2matrix(0.0, 0.0, 0.0) == Matrix(1.0I, 3, 3)
