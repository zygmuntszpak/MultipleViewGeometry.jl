using MultipleViewGeometry, Images
using BenchmarkTools, Compat
using StaticArrays, MAT, Plots

# Load MATLAB matrices that represent a pair of images and that contain
# a set of manually matched corresponding points.
file = matopen("./data/sene-adelaideRMF.mat")
img1 = colorview(RGB, normedview(permutedims(read(file,"img1"),[3,1,2])))
img2 = colorview(RGB, normedview(permutedims(read(file,"img2"),[3,1,2])))
inlier_pts1 = read(file,"inlierPts1")
inlier_pts2 = read(file,"inlierPts2")
close(file)

# Conver the images to Grayscale.
img1g = Gray.(img1)
img2g = Gray.(img2)

plotlyjs()

_, npts1 = size(inlier_pts1)
_, npts2 = size(inlier_pts2)

# Construct a set of corresponding points expressed in homogeneous coordinates.
ℳ = [ Point2DH(vcat(inlier_pts1[:,n],1)) for n in 1:npts1]
ℳʹ =  [ Point2DH(vcat(inlier_pts2[:,n],1)) for n in 1:npts2]

# Estimate the Fundamental matrix from the corresponding points using the
# Direct Linear Transform algorithm.
𝐅₀ = estimate(FundamentalMatrix(),DirectLinearTransform(),  (ℳ, ℳʹ))

# Plot Keypoints, epipole and concomitant epipolar lines in the first view.
p1 = Plots.plot(img1g,grid = false, box = :none, legend = false, size = (455,341))
for n = 1:25:npts1
    m = ℳ[n]
    # Epipolar line in the first image.
    l = 𝑛(𝐅₀'*ℳʹ[n])
    draw!(EpipolarLineGraphic(), l, size(img1), p1)
    # Keypoint.
    Plots.plot!([m[1]],[m[2]], grid = false, box = :none, legend = false,
                    seriestype = :scatter, w = 5, aspect_ratio = :equal)
end
e = epipole(𝐅₀)
Plots.plot!([e[1]],[e[2]], grid = false, box = :none, legend = false,
            seriestype = :scatter, markershape = :diamond, markerstrokecolor = :red,
            markersize = 5, aspect_ratio = :equal)

# Plot Keypoints and concomitant epipolar lines in the second view.
p2 = Plots.plot(img2g,grid = false, box = :none, legend = false)
for n = 1:25:npts2
    mʹ = ℳʹ[n]
    # Epipolar line in the second image.
    l = 𝑛(𝐅₀*ℳ[n])
    draw!(EpipolarLineGraphic(), l, size(img2), p2)
    Plots.plot!([mʹ[1]],[mʹ[2]], grid = false, box = :none, legend = false,
                    seriestype = :scatter, w = 5, aspect_ratio = :equal)
end
eʹ = epipole(𝐅₀')
Plots.plot!([eʹ[1]],[eʹ[2]], grid = false, box = :none, legend = false,
            seriestype = :scatter, markershape = :diamond, markerstrokecolor = :red,
            markersize = 5, aspect_ratio = :equal)


# Display both views simultaneously.
p3 = Plots.plot(p1,p2,layout=(1,2), legend = false)
display(p3)

𝐅₀*e
