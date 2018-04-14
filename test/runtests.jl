using MultipleViewGeometry
using Base.Test


@testset "Operators Tests" begin include("operators_tests.jl") end
@testset "Data Normalization Tests" begin include("data_normalization_tests.jl") end
@testset "Transform Tests" begin include("transform_tests.jl") end
@testset "Fundamental Matrix Construction Tests" begin include("construct_fundamentalmatrix_tests.jl") end
@testset "Estimate Two View Tests" begin include("estimate_twoview_tests.jl") end
@testset "Projection Matrix Construction Tests" begin include("construct_projectionmatrix_tests.jl") end
@testset "Rotations Tests" begin include("rotations_tests.jl") end
@testset "Cost Function Tests" begin include("cost_functions_tests.jl") end
