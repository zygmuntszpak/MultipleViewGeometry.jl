abstract type AbstractVisualProperties end

Base.@kwdef struct  MakieVisualProperties{T₁ <: Real, T₂ <: Number, T₃ <: Number} <: AbstractVisualProperties
    scale::T₁ = 100.0f0
    linewidth::T₂  = 4
    markersize::T₃ = 5
end

Base.@kwdef struct PGFPlotsVisualProperties{T₁ <: Real, T₂ <: Number, T₃ <: Number} <: AbstractVisualProperties
    scale::T₁ = 100.0f0
    linewidth::T₂  = 4
    markersize::T₃ = 5
end
