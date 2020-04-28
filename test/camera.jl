@testset "Camera Model Instantiation" begin
    digital_camera = Camera()
    @unpack model, image_type = digital_camera
    @test typeof(model) <: typeof(Pinhole())
    @test typeof(image_type) <: typeof(AnalogueImage())
end
