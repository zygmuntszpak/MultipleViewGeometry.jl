module ModuleFactorization
using MultipleViewGeometry, MultipleViewGeometry.ModuleTypes
using StaticArrays, LinearAlgebra
export deconstruct, compose
include("deconstruct.jl")
include("compose.jl")
end
