function draw!(camera::CameraModel, scene::AbstractScene)
    optical_center = Point3f0(camera.𝐜)
    image_width = camera.image_width
    image_height = camera.image_height
    f = camera.focal_length
    𝐞₁ = camera.𝐞₁
    𝐞₂ = camera.𝐞₂
    𝐞₃ = camera.𝐞₃
    bottom_right = optical_center + Point3f0((image_width/2)     * 𝐞₁ + (image_height/2)  * 𝐞₂ + f*𝐞₃)
    top_right =  optical_center   + Point3f0((image_width/2)     * 𝐞₁ + (-image_height/2) * 𝐞₂ + f*𝐞₃)
    top_left = optical_center     + Point3f0((-image_width/2)    * 𝐞₁ + (-image_height/2) * 𝐞₂ + f*𝐞₃)
    bottom_left = optical_center  + Point3f0((-image_width/2)    * 𝐞₁ + (image_height/2)  * 𝐞₂ + f*𝐞₃)

    centroid2film = [
        optical_center  => bottom_right;
        optical_center  => top_right;
        optical_center  => top_left;
        optical_center   => bottom_left;
    ]

    film = [
            bottom_right => top_right;
            top_right => top_left;
            top_left => bottom_left;
            bottom_left =>  bottom_right;
            ]

    scale = 20.0f0
    coordinate_system = [
        optical_center => optical_center + Point3f0(scale*𝐞₁);
        optical_center => optical_center + Point3f0(scale*𝐞₂);
        optical_center => optical_center + Point3f0(scale*𝐞₃);
    ]
    linesegments!(scene, centroid2film, color = :black, linewidth = 2)
    linesegments!(scene, film, color = :black, linewidth = 2)
    linesegments!(scene, coordinate_system, color = [:red, :green, :blue ], linewidth = 2)



end
