import AllenNeuropixels as AN
using DataFrames
using Test

@testset "AllenNeuropixels" begin
    cs = @test_nowarn AN.getprobes()
    @test cs isa DataFrame
end
