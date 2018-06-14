abstract type GraphicEntity end

type EpipolarLineGraphic <: GraphicEntity
end


function draw!(g::EpipolarLineGraphic, l::AbstractVector, dim::Tuple{<:Number,<:Number}, p::RecipesBase.AbstractPlot{<:RecipesBase.AbstractBackend})

    top = intersection(l, [0 ; -1  ; 0])
    bottom = intersection(l, [0 ; -1 ; dim[1]])
    left = intersection(l, [-1 ; 0 ; 0])
    right = intersection(l, [-1 ; 0 ; dim[2]])

    x = Float64[]
    y = Float64[]

    if is_inbounds(top,dim)
        push!(x, top[1])
        push!(y, top[2])
    end

    if is_inbounds(bottom,dim)
        push!(x, bottom[1])
        push!(y, bottom[2])
    end

    if is_inbounds(left,dim)
        push!(x, left[1])
        push!(y, left[2])
    end

    if is_inbounds(right,dim)
        push!(x, right[1])
        push!(y, right[2])
    end

    Plots.plot!(x,y,w=3)
end

function intersection(l1::AbstractArray, l2::AbstractArray)
    l = 𝑛(cross(l1,l2))
    l[1:2]
end

function is_inbounds(pt::AbstractVector, dim::Tuple{<:Number,<:Number})
    nrow, ncol = dim
    pt[1] >= -1.5 && pt[1] < ncol+1.5 && pt[2] >= -1.5 && pt[2] <= nrow + 1.5
end
