struct AquireImage <: AbstractContext end



# TODO reconcile this with the ability to change reference coordinate systems (i.e. making a particular camera the coordinate system)
# TODO perhaps allow specifciation of what coordinate system one is considering for the image aquisition
function (aquire::AquireImage)(world::AbstractWorld, camera::AbstractCamera)
    @unpack points = world
    image_points = project(Projection(camera), points)
    return image_points
end
