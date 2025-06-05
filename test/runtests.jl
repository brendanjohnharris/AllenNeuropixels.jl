import AllenNeuropixels as AN
using DataFrames
using Test

@testset "AllenNeuropixels" begin
    cs = AN.getprobes()
    @test cs isa DataFrame
end
