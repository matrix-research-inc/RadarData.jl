#=
    License ID: SEAL_B
=#
export synth_chirp, synth_chirp!, range_compress, range_compress!

"""
    synth_chirp(sample_rate_hz, bandwidth_hz, duration_ns, carrier_hz = 0, start_phase_normalized = 0)

Synthesize chirp given sample rate in Hz, bandwidth in Hz, and duration of the waveform in nanoseconds.
`start_phase_normalized` is a phase offset in cycles, in the range [0,1).
Input a negative `bandwidth_hz` to return a down-chirp.
"""
function synth_chirp(
    sample_rate_hz::Float64,
    bandwidth_hz::Float64,
    duration_ns::Float64,
    carrier_hz::Float64 = 0.0,
    start_phase_normalized::Float64 = 0.0
)
    num_samples = ceil(Int64, sample_rate_hz * duration_ns * 1e-9)
    chirp = zeros(ComplexF64, num_samples)
    synth_chirp!(
        chirp,
        sample_rate_hz,
        bandwidth_hz,
        duration_ns,
        carrier_hz,
        start_phase_normalized
    )
    return chirp
end

"""
    synth_chirp!(chirp, sample_rate_hz, bandwidth_hz, duration_ns, carrier_hz = 0, start_phase_normalized = 0)

Synthesize chirp given sample rate in Hz, bandwidth in Hz, and duration of the waveform in nanoseconds.
`start_phase_normalized` is a phase offset in cycles, in the range [0,1).
Input a negative `bandwidth_hz` to return a down-chirp.
"""
function synth_chirp!(
    chirp::Vector{Complex{T}},
    sample_rate_hz::Float64,
    bandwidth_hz::Float64,
    duration_ns::Float64,
    carrier_hz::Float64 = 0.0,
    start_phase_normalized::Float64 = 0.0
) where {T<:AbstractFloat}
    num_samples = ceil(Int64, sample_rate_hz * duration_ns * 1e-9)
    dt = 1 / sample_rate_hz
    start_freq_hz = carrier_hz - bandwidth_hz/2
    half_chirp_rate = bandwidth_hz / (num_samples * dt * 2)
    chirp .= 0
    for idx in 1:num_samples
        time_s = (idx - 1) * dt
        chirp[idx] = Complex{T}(cis(2 *
                       pi *
                       time_s *
                       (half_chirp_rate * time_s + start_freq_hz) +
                       start_phase_normalized))
    end
    return nothing
end

"""
    range_compress(iq_data, ref_waveform_fft; downsample_rate = 1)

Apply FFT-based convolution to columns of `iq_data` with FFT of a reference waveform. 
"""
function range_compress(
    iq_data::CM,
    ref_waveform_fft::CV,
    downsample_rate::Int = 1
) where {
    T<:AbstractFloat,
    CV<:AbstractVector{Complex{T}},
    CM<:AbstractMatrix{Complex{T}}
}
    if size(iq_data, 1) > length(ref_waveform_fft)
        error("size(iq_data,1) must be less than length(ref_waveform_fft).")
    end
    # Pad the data out
    iq_pad = zeros(Complex{T}, length(ref_waveform_fft), size(iq_data, 2))
    start_pad = fld(length(ref_waveform_fft) - size(iq_data, 1), 2)
    iq_pad[start_pad .+ 1 .+ (1:size(iq_data, 1)), :] = iq_data
    # Transform iq_data to frequency domain 
    iq_fft = fft(iq_pad, 1)
    # Apply conjugate waveform 
    iq_fft .*= conj.(ref_waveform_fft)
    # IFFT back to time domain and optionally downsample
    ifft!(iq_fft, 1)
    if downsample_rate > 1.0
        # Just decimate - may need to revisit?
        mf_data =
            iq_fft[start_pad .+ 1 .+ (1:downsample_rate:size(iq_data, 1)), :]
    else
        mf_data = iq_fft[start_pad .+ 1 .+ (1:size(iq_data, 1)), :]
    end
    return mf_data
end

"""
    range_compress!(iq_data, ref_waveform_fft)

Apply FFT-based convolution to columns of `iq_data` with FFT of a reference waveform. 
This version operates in-place, and so does not offer down-sampling or appropriate padding and trucation.
"""
function range_compress!(
    iq_data::CM,
    ref_waveform_fft::CV
) where {
    T<:AbstractFloat,
    CV<:AbstractVector{Complex{T}},
    CM<:AbstractMatrix{Complex{T}}
}
    # Transform iq_data to frequency domain 
    fft!(iq_data, 1)
    # Apply conjugate waveform 
    iq_data .*= conj.(ref_waveform_fft)
    # IFFT back to time domain
    ifft!(iq_data, 1)
    return nothing
end
