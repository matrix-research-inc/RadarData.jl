@testset "FFT Shortcuts" begin
    @test fft_spacing(3) == -1:1
    @test fft_spacing(4) == -2:1
end
