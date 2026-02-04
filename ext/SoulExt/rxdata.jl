export synth_chirp_fft, VPH

"""
    synth_chirp_fft(tx_msg, rx_msg, sample_rate_hz = 2e9)

Generate the FFT of a chirp waveform defined in the input Tx and Rx messages.
Provides a default padding value for FFT-based filtering.
"""
function synth_chirp_fft(
    tx_msg::DmaTxMsg,
    rx_msg::DmaRxMsg,
    sample_rate_hz = 2e9
)
    # Infer dimensions of receive data and choose padding 
    tx_samples = ceil(Int64, tx_msg.duration_ns * 1e-9 * sample_rate_hz)
    num_samples = ceil(Int64, rx_msg.duration_ns * 1e-9 * sample_rate_hz)
    fft_size = nextprod([2, 3, 5, 7, 11], num_samples + max(tx_samples, 200))

    adjust_center_freq_hz =
        tx_msg.start_freq_baseband_hz + tx_msg.bandwidth_hz / 2
    chirp_fft = zeros(ComplexF32, fft_size)
    synth_chirp!(
        chirp_fft,
        sample_rate_hz,
        tx_msg.bandwidth_hz,
        tx_msg.duration_ns,
        adjust_center_freq_hz,
        tx_msg.starting_phase_normalized
    )
    fft!(chirp_fft)
    return chirp_fft
end

"""
    VPH(rx_data, tx_msg, rx_msg; carrier_hz = 3e9, sample_rate_hz = 2e9, calibration_offset_ns = 0, filter_fft = synth_chirp_fft())

Construct a video phase history from SoulHardware's RxData type, and metadata provided by 
the transmit and receive messages sent to the radar hardware.
NOTE: The `ref_time_s` field will be populated with the number of seconds since FPGA power-on, 
not the usual time since Unix Epoch which is normally provided in a VPH.

# Arguments:
- `center_freq_hz`: Center (or carrier) frequency in Hz 
- `sample_rate_hz`: Complex sample rate of ADC, in Hz 
- `calibration_offset_ns`: Timing in ns it takes for signals to travel through RF chain to antenna and back
- `filter_fft`: Fourier transform of the desired waveform for filtering the IQ data in `rx_data`
"""
function RadarData.VPH(
    rx_data::RxData,
    tx_msg::DmaTxMsg,
    rx_msg::DmaRxMsg;
    center_freq_hz::Float64 = 3e9,
    sample_rate_hz::Float64 = 2e9,
    calibration_offset_ns::Float64 = 0.0,
    filter_fft::Vector{ComplexF32} = synth_chirp_fft(
        tx_msg,
        rx_msg,
        sample_rate_hz
    ),
    presum_factor::Int64 = 1
)::RadarData.VPH{Float32}
    # Infer data size 
    num_pulses = fld(length(rx_data.frames), presum_factor)
    num_samples = Int64(length(rx_data.frames[1]))

    freq_list_hz =
        center_freq_hz .+
        collect(fftshift(fftfreq(num_samples, sample_rate_hz)))
    # Delay to first range bin
    two_way_delay_s =
        (rx_msg.start_time_ns - tx_msg.start_time_ns - calibration_offset_ns) *
        1e-9
    # Delay to center range bin
    center_delay = fld(num_samples, 2) / sample_rate_hz
    ref_range_m = c0 * (two_way_delay_s + center_delay) .* ones(num_pulses)
    slow_time_s = collect(range(
        0.0,
        length = num_pulses,
        step = rx_msg.pri_spacing_ns * 1e-9 * presum_factor
    ))

    # This is time since FPGA on, in clock cycles (2 ns)
    # The standard VPH field is seconds since Unix Epoch, so this
    # will be modified downstream if necessary
    ref_time_s = Float64(rx_data.header.system_time)

    # Optionally pre-sum pulses in batches
    reshape_rx = zeros(ComplexF32, num_samples, num_pulses)
    if presum_factor > 1
        for cidx in axes(reshape_rx, 2)
            idxs = ((cidx - 1) * presum_factor .+
             (1:presum_factor))::UnitRange{Int64}
            reshape_rx[:, cidx] = sum(rx_data.frames[idxs])
        end
        reshape_rx ./= presum_factor
    else
        reshape_rx[:, :] = reshape(rx_data.buffer, num_samples, num_pulses)
    end

    # Apply range compression
    mf_data = range_compress(reshape_rx, filter_fft)

    # Placeholders
    tx_pos_enu = zeros(3, 1)
    rx_pos_enu = zeros(3, 1)
    srp_enu = zeros(3, 1)
    srp_lla = SVector{3,Float64}(0, 0, 0)

    return RadarData.VPH(
        mf_data,
        freq_list_hz,
        ref_range_m,
        slow_time_s,
        ref_time_s,
        tx_pos_enu,
        rx_pos_enu,
        srp_enu,
        srp_lla,
        -1.0,
        Float64(c0),
        true
    )
end

# """
#     load_rxdata(file_name)

# Read an RxData object, along with system parameters from a BEVE file.
# """
# function load_rxdata(file_name::String)
#     reader = BeveDeserializer(open(file_name, "r"))
#     rx_beve = parse_value(reader)

#     headers = rx_beve["rx_data"][1][1]["headers"]
#     frames = rx_beve["rx_data"][1][1]["frames"]

#     sample_rate_hz =
#         rx_beve["meta"]["system"]["clock_rates"]["complex_sample_rate_hz"]

#     dma_msgs = rx_beve["meta"]["dma_msgs"]

#     # TODO: There *should* only be one of each but who knows.
#     tx_msg_ind = findfirst(msg -> msg.index == 0, dma_msgs)
#     rx_msg_ind = findfirst(msg -> msg.index == 1, dma_msgs)

#     tx_msg = dma_msgs[tx_msg_ind].value
#     rx_msg = dma_msgs[rx_msg_ind].value

#     return rx_data, tx_msg, rx_msg, sample_rate_hz
# end