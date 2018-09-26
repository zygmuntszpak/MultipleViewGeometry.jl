module ModuleConstraints
using MultipleViewGeometry.ModuleTypes, MultipleViewGeometry.ModuleMathAliases, MultipleViewGeometry.ModuleOperators
using StaticArrays, LinearAlgebra
export satisfy, EpipolarConstraint, Constraint

include("constraints.jl")
end
