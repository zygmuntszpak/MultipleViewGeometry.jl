using MultipleViewGeometry
using LinearAlgebra
using Test
using StaticArrays
using GeometryTypes
using Parameters

@testset "MultipleViewGeometry.jl" begin
    @testset "Struct Instantiation"   begin
            include("camera.jl")
            include("planar_scene.jl")
    end
    @testset "Intrinsics"   begin
            include("intrinsics_coordinate_systems.jl")
    end
    @testset "Fundamental Matrix"   begin
            include("fundamental_matrix.jl")
    end
    @testset "Homography Matrix"   begin
            include("homography_matrix.jl")
    end
end
