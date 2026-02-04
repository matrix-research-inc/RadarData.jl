#=
    License ID: SEAL_B
=#
using Test
using RadarData
using DSP: conv
using FFTW

# Synthesize a chirp 
sample_rate_hz = 2.5e9
bandwidth_hz = 100e6
duration_ns = 100.0
chirp = synth_chirp(sample_rate_hz, bandwidth_hz, duration_ns)

# using Plots
# plot(
#     fft_spacing(length(chirp)) * sample_rate_hz / length(chirp) * 1e-6,
#     abs.(fftshift(fft(chirp)))
# )
# title!("Spectrum of Chirp")
# xlabel!("Frequency, MHz")

# Generate some sample data

sample_locs = zeros(ComplexF64, 1000, 5)
sample_locs[500, 1] = 1
sample_locs[501, 2] = 1
sample_locs[502, 3] = 1
sample_locs[503, 4] = 1
sample_locs[504, 5] = 1

iq_data = conv(sample_locs, chirp)

# Test range compression
fft_size = nextprod([2, 3, 5, 7], size(iq_data, 1) + 100)
ref_waveform_fft = zeros(ComplexF64, fft_size)
ref_waveform_fft[1:length(chirp)] = chirp
fft!(ref_waveform_fft)

compressed = range_compress(iq_data, ref_waveform_fft)

# Peaks should be in the same location
@test argmax(abs.(compressed[:, 1])) == 500
@test argmax(abs.(compressed[:, 2])) == 501
@test argmax(abs.(compressed[:, 3])) == 502
@test argmax(abs.(compressed[:, 4])) == 503
@test argmax(abs.(compressed[:, 5])) == 504

# in-place, when types are appropriate
ref_waveform_fft = zeros(ComplexF64, size(iq_data, 1))
ref_waveform_fft[1:length(chirp)] = chirp
fft!(ref_waveform_fft)
range_compress!(iq_data, ref_waveform_fft)

@test argmax(abs.(iq_data[:, 1])) == 500
@test argmax(abs.(iq_data[:, 2])) == 501
@test argmax(abs.(iq_data[:, 3])) == 502
@test argmax(abs.(iq_data[:, 4])) == 503
@test argmax(abs.(iq_data[:, 5])) == 504
