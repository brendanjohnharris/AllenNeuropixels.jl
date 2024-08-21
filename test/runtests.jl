@testset "AllenNeuropixels" begin
    import AllenNeuropixels as AN
    using DataFrames
    cs = @test_nowarn AN.getprobes()
    @test cs isa DataFrame
end
