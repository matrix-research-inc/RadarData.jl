using RadarData
using Test
using Aqua

@testset "RadarData.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(RadarData)
    end

    include("fft_shortcuts.jl")
end
