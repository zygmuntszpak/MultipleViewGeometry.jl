using MultipleViewGeometry
using LinearAlgebra
using Test
using StaticArrays
using GeometryBasics
using Parameters
using Random
using FiniteDiff

include("util.jl")

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
    @testset "Essential Matrix"   begin
            include("essential_matrix.jl")
    end
    @testset "Homography Matrix"   begin
            include("homography_matrix.jl")
    end

    @testset "Camera Calibration"   begin
            include("sole_camera_rig_calibration.jl")
    end

    @testset "Pose and Extrinsics"   begin
            include("pose_extrinsics.jl")
    end

    @testset "Coordinate Transformations"   begin
            include("world_system_transformation.jl")
    end

    @testset "Triangulation"   begin
            include("triangulate.jl")
    end
end
