module ModuleCostFunction
using MultipleViewGeometry.ModuleTypes, MultipleViewGeometry.ModuleOperators, MultipleViewGeometry.ModuleMathAliases
using MultipleViewGeometry.ModuleCarriers
using MultipleViewGeometry.ModuleTransform
using StaticArrays, LinearAlgebra
export cost, X, H, covariance_matrix, covariance_matrix_debug, covariance_matrix_normalised
include("cost_functions.jl")
end
