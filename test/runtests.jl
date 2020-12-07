using PBDS
using Test

@testset "PBDS.jl" begin
    @test dim(S1) == 1
    @test embdim(S1) == 2
end
