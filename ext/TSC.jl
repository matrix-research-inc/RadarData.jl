"""An extension to RadarData for loading Technology Service Corporation data files into VPH format."""
module TSC

using MAT: matread
using FFTW: fft!
using RadarData: c0, VPH, range_compress
using Geodesy: ENU, LLA, ECEF, wgs84

const leapsec_dates_posix = [
    362793600,
    394329600,
    425865600,
    489024000,
    567993600,
    631152000,
    662688000,
    709948800,
    741484800,
    773020800,
    820454400,
    867715200,
    915148800,
    1136073600,
    1230768000,
    1341100800,
    1435708800,
    1483228800
]

"""
    gps_to_utc(gps_seconds, gps_week = 0)

Convert GPS time to UTC, accounting for leap seconds.
Returns a POSIX time, seconds since Unix epoch.
"""
function gps_to_utc(gps_seconds::Float64, gps_week::Int = 0)
    diff_gps_minus_unix_epoch = 315964800
    gps_seconds += gps_week * 604800
    base_utc_time = gps_seconds + diff_gps_minus_unix_epoch
    for date in leapsec_dates_posix
        if base_utc_time > date
            base_utc_time += 1
        end
    end
    return base_utc_time
end

"""
    load_tsc(file_name, gps_week = 0; sub_inds)

Load data from Technology Service Corporation (TSC) in .mat format.
Applies range compression.
First data drop did not include GPS week, so that must be passed manually for now.   
"""
function load_tsc(
    file_name::String,
    gps_week::Int = 0;
    sub_inds::AV = [0]
) where {AV<:AbstractVector{Int}}
    mat_file = matread(file_name)
    c_eff_ms = Float64(c0)

    # Load IQ data
    iq_data = mat_file["iq"]
    if any(sub_inds .< 1)
        sub_inds = 1:size(iq_data, 2)
    end
    iq_data = iq_data[:, sub_inds]
    num_samples, num_pulses = size(iq_data)

    # Extract metadata 
    waveform = mat_file["systemParameters"]["refWaveform"]
    center_freq_hz = mat_file["systemParameters"]["txFrequency_khz"] * 1e3
    # Real-valued sampling rate, so exactly twice bandwidth
    bandwidth_hz = mat_file["systemParameters"]["txBandwidth_hz"]
    sampling_rate_hz = mat_file["systemParameters"]["samplingRate_hz"]
    range_step_m = c_eff_ms / bandwidth_hz / 2
    range_first_m =
        c_eff_ms * 1e-9 / 2 * mat_file["systemParameters"]["rxDelay_ns"]
    # one-way range
    ref_range_m = range_first_m + num_samples / 2 * range_step_m

    # Convert timestamps from relative GPS time to UTC time
    if gps_week == 0
        @warn "No GPS reference time given. UTC time may be incorrect!"
    end
    time_stamps = gps_to_utc.(
        vec(mat_file["pulseTimes"]["processedTimestamp"])[sub_inds],
        gps_week
    )

    # Tx path and SRP estimation
    center_pulse = fld(num_pulses, 2) + 1
    tx_pos_lla = mat_file["ins"]["lla"][:, sub_inds]
    # Antenna mounted perpendicular to body 
    # Approximately 5 or 4.3 degree depression
    look_dir_mult = mat_file["systemParameters"]["lookDir"] == "right" ? 1 : -1
    pra = mat_file["ins"]["PRA"][:, sub_inds] # pitch-roll-azimuth
    # Antenna angles from North
    thetas = pi / 180 * (pra[3, :] .+ look_dir_mult .* 90)
    # Estimate distance in ground plane to SRP
    ground_ranges = sqrt.(ref_range_m .^ 2 .- tx_pos_lla[3, :] .^ 2)
    # North and East line of sight vector
    los_vecs_enu = hcat(sin.(thetas), cos.(thetas), zeros(size(thetas)))
    # Estimate SRP as point where LOS intersects ground plane
    srp_rels = transpose(ground_ranges .* los_vecs_enu)
    srp_rels[3, :] = -tx_pos_lla[3, :]
    # Scene reference point for collect
    srp_ref_ecef = ECEF(
        ENU(
            srp_rels[1, center_pulse],
            srp_rels[2, center_pulse],
            srp_rels[3, center_pulse]
        ),
        ECEF(
            LLA(
                tx_pos_lla[1, center_pulse],
                tx_pos_lla[2, center_pulse],
                tx_pos_lla[3, center_pulse]
            ),
            wgs84
        ),
        wgs84
    )
    srp_lla = LLA(srp_ref_ecef, wgs84)
    srp_enu = zeros(3, num_pulses)
    tx_pos_enu = zeros(3, num_pulses)
    for pidx in 1:num_pulses
        local tx_pos_ecef = ECEF(
            LLA(tx_pos_lla[1, pidx], tx_pos_lla[2, pidx], tx_pos_lla[3, pidx]),
            wgs84
        )
        # Compute ENU coordinates for flight path
        tx_pos_enu[:, pidx] = ENU(tx_pos_ecef, srp_ref_ecef, wgs84)
        local srp_ecef = ECEF(
            ENU(srp_rels[1, pidx], srp_rels[2, pidx], srp_rels[3, pidx]),
            tx_pos_ecef,
            wgs84
        )
        # Relative SRP (for stripmap)
        srp_enu[:, pidx] = ENU(srp_ecef, srp_ref_ecef, wgs84)
    end

    # Apply matched filter to data 
    fft_size = nextprod([2, 3, 5, 7, 9], num_samples + 200)
    ref_waveform_fft = zeros(ComplexF64, fft_size)
    ref_waveform_fft[1:length(waveform)] = waveform ./ 840
    fft!(ref_waveform_fft)
    mf_data = range_compress(iq_data, ref_waveform_fft)

    freq_list_hz = collect(range(
        center_freq_hz - bandwidth_hz / 2,
        center_freq_hz + bandwidth_hz / 2,
        size(mf_data, 1)
    ))

    kernel_sign = -1.0
    range_domain = true
    return VPH(
        mf_data,
        freq_list_hz,
        fill(2.0 * ref_range_m, num_pulses),
        time_stamps .- time_stamps[1],
        time_stamps[1],
        tx_pos_enu,
        tx_pos_enu,
        srp_enu,
        [srp_lla.lat, srp_lla.lon, srp_lla.alt],
        kernel_sign,
        c_eff_ms,
        range_domain
    )
end

end
