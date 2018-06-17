abstract type GraphicEntity end

type EpipolarLineGraphic <: GraphicEntity
end

type LineSegment3D <: GraphicEntity
end

type PlaneSegment3D <: GraphicEntity
end

type Camera3D <: GraphicEntity
end

type WorldCoordinateSystem3D <: GraphicEntity
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


function draw!(g::LineSegment3D, 𝐨::AbstractArray, 𝐩::AbstractArray, col::Symbol, p::RecipesBase.AbstractPlot{<:RecipesBase.AbstractBackend})
    x = [𝐨; 𝐩][:,1]
    y = [𝐨; 𝐩][:,2]
    z = [𝐨; 𝐩][:,3]
    Plots.path3d!(x,y,z, w = 2,grid = false, box = :none, legend = false, linecolor = col)
end

function draw!(g::LineSegment3D, 𝐨::AbstractVector, 𝐩::AbstractVector, col::Symbol, p::RecipesBase.AbstractPlot{<:RecipesBase.AbstractBackend})
    draw!(LineSegment3D(), 𝐨', 𝐩', col, p)
end

function draw!(g::PlaneSegment3D,  𝐩₁::AbstractArray, 𝐩₂::AbstractArray, 𝐩₃::AbstractArray, 𝐩₄::AbstractArray, col::Symbol, p::RecipesBase.AbstractPlot{<:RecipesBase.AbstractBackend})
    draw!(LineSegment3D(), 𝐩₁, 𝐩₂, col, p)
    draw!(LineSegment3D(), 𝐩₂, 𝐩₃, col, p)
    draw!(LineSegment3D(), 𝐩₃, 𝐩₄, col, p)
    draw!(LineSegment3D(), 𝐩₄, 𝐩₁, col, p)
end

function draw!(g::WorldCoordinateSystem3D, scale,  p::RecipesBase.AbstractPlot{<:RecipesBase.AbstractBackend})
    𝐞₁ = [1,  0,  0]
    𝐞₂ = [0,  1,   0]
    𝐞₃ = [0,  0,   1]
    𝐨  = [0,  0,  0]

    # Draw the world coordinate axes.
    draw!(LineSegment3D(), 𝐨, 𝐨 + scale*𝐞₁, :red, p)
    draw!(LineSegment3D(), 𝐨, 𝐨 + scale*𝐞₂, :green, p)
    draw!(LineSegment3D(), 𝐨, 𝐨 + scale*𝐞₃, :blue, p)
end

function draw!(g::Camera3D, 𝐊::AbstractArray,  𝐑::AbstractArray, 𝐭::AbstractArray, scale,  p::RecipesBase.AbstractPlot{<:RecipesBase.AbstractBackend})

    # Origin of the world coordinate system.
    𝐞₁ = [1,  0,  0]
    𝐞₂ = [0,  1,   0]
    𝐞₃ = [0,  0,   1]
    𝐨  = [0,  0,  0]

    # Initial camera imaging plane.
    𝐩₁ =  [-125,  125,  -50]
    𝐩₂ =  [125,  125,  -50]
    𝐩₃ =  [125, -125, -50]
    𝐩₄ =  [-125,  -125, -50]

    # Initial camera center.
    𝐜 = [0.0, 0.0, 0.0]
    𝐜 = 𝐑*𝐜 + 𝐭
    Plots.plot!([𝐜[1]],[𝐜[2]],[𝐜[3]],seriestype = :scatter, ms=1, grid = false, box = :none, legend = false, markercolor=:Red)

    draw!(PlaneSegment3D(), 𝐑*𝐩₁ + 𝐭, 𝐑*𝐩₂ + 𝐭, 𝐑*𝐩₃ + 𝐭, 𝐑*𝐩₄ + 𝐭, :black, p)

    # Connect camera center with corners of plane segment.
    draw!(LineSegment3D(), 𝐜, 𝐑*𝐩₁ + 𝐭, :black, p)
    draw!(LineSegment3D(), 𝐜, 𝐑*𝐩₂ + 𝐭, :black, p)
    draw!(LineSegment3D(), 𝐜, 𝐑*𝐩₃ + 𝐭, :black, p)
    draw!(LineSegment3D(), 𝐜, 𝐑*𝐩₄ + 𝐭, :black, p)

    # Draw camera coordinate axes for the first camera.
    draw!(LineSegment3D(), 𝐜, (𝐑*scale*𝐞₁ + 𝐭), :red, p)
    draw!(LineSegment3D(), 𝐜, (𝐑*scale*𝐞₂ + 𝐭), :green, p)
    draw!(LineSegment3D(), 𝐜, (𝐑*scale*𝐞₃ + 𝐭), :blue, p)


end
