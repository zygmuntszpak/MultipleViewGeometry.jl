using MultipleViewGeometry
using Base.Test


@testset "Operators Tests" begin include("operators_tests.jl") end
@testset "Data Normalization Tests" begin include("data_normalization_tests.jl") end
