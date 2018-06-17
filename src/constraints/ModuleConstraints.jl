module ModuleConstraints
using MultipleViewGeometry.ModuleTypes, MultipleViewGeometry.ModuleMathAliases, MultipleViewGeometry.ModuleOperators
using StaticArrays
export satisfy, EpipolarConstraint, Constraint

include("constraints.jl")
end
