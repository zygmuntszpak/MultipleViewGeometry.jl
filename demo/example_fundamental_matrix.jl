using MultipleViewGeometry, Images, ImageFeatures, ImageView
using BenchmarkTools, Compat
using StaticArrays, MAT, Plots, PyPlot

file = matopen("./data/sene-adelaideRMF.mat")
# img1 = read(file,"img1")
# img2 = read(file,"img2")
img1 = colorview(RGB, normedview(permutedims(read(file,"img1"),[3,1,2])))
img2 = colorview(RGB, normedview(permutedims(read(file,"img2"),[3,1,2])))
inlier_pts1 = read(file,"inlierPts1")
inlier_pts2 = read(file,"inlierPts2")
close(file)

img1g = Gray.(img1)
img2g = Gray.(img2)

x = 1:10; y = rand(10); # These are the plotting data
plot(x,y)
#pyplot()
plotly()

p = Plots.plot(img1g)
Plots.plot!(inlier_pts1[1,:],inlier_pts1[2,:],seriestype=:scatter,w=5)

p2 = Plots.plot(img2g)
Plots.plot!(inlier_pts2[1,:],inlier_pts2[2,:],seriestype=:scatter,w=5)

#Plots.plot(img1g)
#Plots.plot!(x->200sin(.05x)+300, 0, 700, seriestype=:scatter,w=5)

#img = reinterpret(N0f8, img1)
# test = permutedims(img1,[3,1,2])
# img_l = colorview(RGB, normedview(test))
# img_r = colorview(RGB, normedview(test))
# img1 = Gray.(img1)
#
# img = colorview(Gray,img1)
#
# colorview(RGB,img1)
#
# img = reinterpret(N0f8, img1)
# test = colorview(RGB,img,size(img))
#
# test = channelview(img)
#
# Gray.(test)
#
# test = permutedims(img,[3,1,2])
# Gray.(channelview(test))
# Gray.(test)
#
# colorview(RGB, normedview(test))


# using PyPlot
# x = linspace(0,2*pi,1000); y = sin.(3 * x + 4 * cos.(2 * x));
# plot(x, y, color="red", linewidth=2.0, linestyle="--")
# title("A sinusoidally modulated sinusoid")
