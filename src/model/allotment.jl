abstract type AbstractAllotment end

struct IntervalAllotment{T₁ <: AbstractRange} <: AbstractAllotment
    interval::T₁
end

# function get_interval(allotment::IntervalAllotment)
#     allotment.interval
# end
