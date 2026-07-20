export VPH,
       get_freq_list,
       spatial_freqs,
       bandwidth,
       center_freq,
       range_res,
       velocity_res,
       range_axis,
       fast_time_axis,
       velocity_axis,
       time_axis,
       doppler_axis,
       cross_range_extent,
       bistatic_angles,
       bisectors,
       vph_to_rs!,
       vph_to_rs,
       rs_to_vph!,
       rs_to_vph,
       freq_truncate!,
       time_truncate!,
       freq_pad!,
       range_pad!,
       extract_cpi,
       pri,
       azimuth_rad,
       elevation_rad,
       horizon_coords,
       range_dir_enu,
       cross_range_dir_enu,
       bp_layover_projection

"""
    Video Phase History data structure.

Each column represents a single pulse of returned data, which has been range-compressed.
The columns may represent either frequency or fast-time data, according to the value of `range_domain`.
Data may be either ComplexF32 or ComplexF64, but position and timing metadata will be retained at Float64 precision.

$(TYPEDFIELDS)
"""
mutable struct VPH{T <: AbstractFloat}
    "Complex-valued data matrix, frequency (or fast-time / range) by slow-time"
    data::Matrix{Complex{T}}
    "Frequency list, in Hz"
    freq_list_hz::Vector{Float64}
    "Two-way range to the center range bin for each pulse"
    ref_range_m::Vector{Float64}
    "Per-pulse collection time, in seconds relative to `ref_time_s`"
    slow_time_s::Vector{Float64}
    "Reference time, in posix time (seconds since Unix Epoch)"
    ref_time_s::Float64
    "Position of transmitter in East-North-Up (ENU), meters"
    tx_pos_enu::Matrix{Float64}
    "Position of receiver (ENU, meters)"
    rx_pos_enu::Matrix{Float64}
    "Position of (possibly shifting) scene reference point (ENU, meters)"
    srp_enu::Matrix{Float64}
    "Fixed scene reference position for ENU frame in latitude (deg), longitude (deg), altitude (m)"
    srp_lla::SVector{3, Float64}
    "FFT sign convention for transformation from frequency to fast-time domain, default is -1"
    kernel_sign::Float64
    "Effective speed of light (m/s)"
    c_eff_ms::Float64
    "Whether the first dimension represents frequency or range (fast-time)"
    range_domain::Bool
end

"""
	VPH(data, freq_list_hz, ref_range_m, slow_time_s, tx_pos_enu; kw_args)

Convenience constructor for VPH type.

Keyword Arguments:
- `ref_time_s = 0`
- `rx_pos_enu = tx_pos_enu`
- `srp_enu = zeros(3,1)`
- `srp_lla = [39.729866; -84.077144; 0]`
- `kernel_sign = -1`
- `c_eff_ms = c0`
- `range_domain = false`
"""
VPH(
data::AbstractMatrix{Complex{T}},
freq_list_hz::AbstractVector,
ref_range_m::AbstractVector,
slow_time_s::AbstractVector,
tx_pos_enu::AbstractMatrix;
ref_time_s::Real = 0.0,
rx_pos_enu::AbstractMatrix = tx_pos_enu,
srp_enu::AbstractMatrix = zeros(3, 1),
srp_lla::AbstractVector = [39.729866; -84.077144; 0.0],
kernel_sign::Real = -1.0,
c_eff_ms::Real = Float64(c0),
range_domain::Bool = false
) where {T <: AbstractFloat} = VPH(
    convert(Matrix{Complex{T}}, data),
    convert(Vector{Float64}, freq_list_hz),
    convert(Vector{Float64}, ref_range_m),
    convert(Vector{Float64}, slow_time_s),
    Float64(ref_time_s),
    convert(Matrix{Float64}, tx_pos_enu),
    convert(Matrix{Float64}, rx_pos_enu),
    convert(Matrix{Float64}, srp_enu),
    convert(SVector{3, Float64}, srp_lla),
    Float64(kernel_sign),
    Float64(c_eff_ms),
    range_domain
)

"""
    VPH(n_freqs, n_pulses)

Convenience constructor for generating a quick test VPH of a given size.
"""
VPH(n_freqs::Int, n_pulses::Int) = VPH(
    zeros(ComplexF64, n_freqs, n_pulses),
    zeros(n_freqs),
    zeros(n_pulses),
    Float64.(collect(1:n_pulses)),
    zeros(3, n_pulses)
)

function Base.convert(::Type{VPH{T2}}, vph::VPH{T1}) where {T1, T2 <: AbstractFloat}
    VPH(
        convert(Matrix{Complex{T2}}, vph.data),
        vph.freq_list_hz,
        vph.ref_range_m,
        vph.slow_time_s,
        vph.ref_time_s,
        vph.tx_pos_enu,
        vph.rx_pos_enu,
        vph.srp_enu,
        vph.srp_lla,
        vph.kernel_sign,
        vph.c_eff_ms,
        vph.range_domain
    )
end

Base.size(vph::VPH) = size(vph.data)
Base.size(vph::VPH, dims::Int) = size(vph.data, dims)
Base.getindex(vph::VPH, idx) = vph.data[idx]
Base.eltype(vph::VPH) = eltype(vph.data)

"""
    Approximate pulse repetition interval of VPH.
    Assumes constant pulse spacing.
"""
pri(vph::VPH) = mean(diff(vph.slow_time_s))

"""
	get_freq_list(center_freq, bandwidth, num_freqs)
"""
get_freq_list(center_freq, bandwidth, num_freqs) = fft_spacing(num_freqs) * bandwidth /
                                                   num_freqs .+ center_freq

"""
	spatial_freqs(freqs)

For frequency `freqs` in Hz, return the spatial frequency `2πf/c`.
"""
spatial_freqs(freqs::AbstractVector{T}, c_eff = c0) where {T <: AbstractFloat} = T(2 * pi /
                                                                                   c_eff) .*
                                                                                 freqs

"""
    spatial_freqs(vph)

For frequency `freqs` in Hz, return the spatial frequency `2πf/c`.
"""
spatial_freqs(vph::VPH) = spatial_freqs(vph.freq_list_hz, vph.c_eff_ms)

"""
	bandwidth(freq_list_hz)

Return the bandwidth represented by equispaced `freq_list_hz`.
"""
bandwidth(freq_list_hz) = (freq_list_hz[2] - freq_list_hz[1]) * length(freq_list_hz)

"""
    bandwidth(vph)

Return the bandwidth spanned by `VPH.freq_list_hz`.
"""
bandwidth(vph::VPH) = bandwidth(vph.freq_list_hz)

"""center_freq(freq_list_hz)"""
center_freq(freq_list_hz) = freq_list_hz[div(length(freq_list_hz), 2) + 1]
"""center_freq(vph)"""
center_freq(vph::VPH) = center_freq(vph.freq_list_hz)

"""
	range_res(bandwidth)

Given `bandwidth` in Hz, compute two-way range resolution in meters.
Note that this value will not change for upsampled range profiles.
"""
range_res(bandwidth, c_eff = c0) = c_eff / bandwidth

"""
	range_res(vph)

Compute two-way range resolution in meters.
"""
range_res(vph::VPH) = range_res(bandwidth(vph), vph.c_eff_ms)

"""
	velocity_res(center_freq)

Given `center_freq` in Hz, compute two-way velocity resolution in meters per pulse.
"""
velocity_res(center_freq, num_pulses, c_eff = c0) = c_eff / center_freq / num_pulses

"""
	velocity_res(vph)

Compute two-way velocity resolution in meters per pulse.
"""
velocity_res(vph::VPH) = velocity_res(center_freq(vph), size(vph, 2), vph.c_eff_ms)

"""
	range_axis(range_step, num_bins)

Compute range axis (meters) of range-Doppler image given the range step and number of range bins.
"""
range_axis(range_step, num_bins) = fft_spacing(num_bins) * range_step

"""
    range_axis(vph)

Compute range axis (two-way meters) of range-Doppler image formed from VPH.
"""
range_axis(vph::VPH) = range_axis(range_res(vph), size(vph, 1))

"""
    fast_time_axis(vph)

Compute the range axis of the VPH, scaled to relative times of arrival, in seconds.
"""
fast_time_axis(vph::VPH) = range_axis(vph) / c0

"""
	velocity_axis(vel_step, num_bins)

Compute velocity axis (m/pulse) of range-Doppler image formed from VPH.
"""
velocity_axis(vel_step, num_bins) = fft_spacing(num_bins) * vel_step

"""
	velocity_axis(vph)

Compute velocity axis (m/pulse) of range-Doppler image formed from VPH.
"""
velocity_axis(vph::VPH) = velocity_axis(velocity_res(vph), size(vph, 2))

"""
    time_axis(vph)

Form an equispaced list of times, spanning duration of VPH collection.
"""
time_axis(vph::VPH) = LinRange(vph.slow_time_s[1], vph.slow_time_s[end], size(vph, 2))

"""
    doppler_axis(vph)

Compute Doppler frequency axis (Hz) for range-Doppler map.
"""
function doppler_axis(vph::VPH)
    prf_limit = 0.5 / pri(vph)
    return LinRange(-prf_limit, prf_limit, size(vph, 2))
end

"""
	cross_range_extent(vph)

Compute unambiguous cross range extent (in meters) around scene reference point.
"""
function cross_range_extent(vph::VPH)
    pulse_az = azimuth_rad(vph)
    delta_az = mean(diff(pulse_az))
    return abs(vph.c_eff_ms / 2 / center_freq(vph) / delta_az)
end

"""
	bistatic_angles(vph)

Compute bistatic angle per-pulse.
"""
function bistatic_angles(vph::VPH)
    inner_prod = dot.(eachcol(vph.tx_pos_enu), eachcol(vph.rx_pos_enu))
    tx_mags = col_norm(vph.tx_pos_enu)
    rx_mags = col_norm(vph.rx_pos_enu)
    scale_prod = clamp.(inner_prod ./ (tx_mags .* rx_mags), -1.0, 1.0)
    return acos.(scale_prod)
end

"""
	bisectors(vph)

Return bistatic bisectors for VPH object.
"""
function bisectors(vph::VPH)
    unit_tx = vph.tx_pos_enu ./ col_norm(vph.tx_pos_enu)
    unit_rx = vph.rx_pos_enu ./ col_norm(vph.rx_pos_enu)
    return unit_tx .+ unit_rx
end

"""
    range_dir_enu(vph)

Return a unit vector pointing towards the down-range direction, in ENU coordinates.
"""
function range_dir_enu(vph::VPH)
    center_pulse = fld(size(vph, 2), 2) + 1
    unit_tx = vph.tx_pos_enu[:, center_pulse] ./
              sqrt(sum(abs2, vph.tx_pos_enu[:, center_pulse]))
    unit_rx = vph.rx_pos_enu[:, center_pulse] ./
              sqrt(sum(abs2, vph.rx_pos_enu[:, center_pulse]))
    bisector = unit_tx .+ unit_rx
    bisector ./= sqrt(sum(abs2, bisector))
    return bisector
end

"""
    cross_range_dir_enu(vph)

Return a unit vector pointing towards the cross-range direction (positive velocity), in ENU coordinates.
"""
function cross_range_dir_enu(vph::VPH)
    center_pulse = fld(size(vph, 2), 2) + 1
    bisec = bisectors(vph)
    # Project out range direction
    range_dir = bisec[:, center_pulse] ./ sqrt(sum(abs2, bisec[:, center_pulse]))
    bisec .-= range_dir .* (range_dir' * bisec)
    # Find remaining axis
    fact = LinearAlgebra.svd(bisec)
    cross_range_dir = fact.U[:, 1]
    # Make sure the vector is pointing towards positive velocity
    rel_pos = cross_range_dir' * bisec
    if mean(diff(vec(rel_pos))) < 0
        cross_range_dir .*= -1
    end
    return cross_range_dir
end

"""
    bp_layover_projection(vph)

Compute approximate layover projection for ENU coordinates in short aperture backprojection imagery.
Only appropriate for approximately monostatic collections.
"""
function bp_layover_projection(vph::VPH)
    cidx = fld(size(vph, 2), 2) + 1
    mono_pos_enu = 0.5 .*
                   (vph.tx_pos_enu[:, cidx .+ (-1:1)] .+ vph.rx_pos_enu[:, cidx .+ (-1:1)])
    mono_velocity_enu = (mono_pos_enu[:, 3] .- mono_pos_enu[:, 1]) ./
                        (vph.slow_time_s[cidx + 1] - vph.slow_time_s[cidx - 1])
    proj_vec = cross(mono_pos_enu[:, 2], mono_velocity_enu)
    return SMatrix{3, 3, Float64}(
        1,
        0,
        0,
        0,
        1,
        0,
        proj_vec[1] / proj_vec[3],
        proj_vec[2] / proj_vec[3],
        0
    )
end

"""Compute azimuth look direction, in radians."""
azimuth_rad(vph::VPH) = azimuth_rad(bisectors(vph))
"""Compute elevation look direction, in radians."""
elevation_rad(vph::VPH) = elevation_rad(bisectors(vph))
"""Compute horizon coordinates (range, azimuth, elevation), in meters and radians."""
horizon_coords(vph::VPH) = cart_to_horizon_unwrap(bisectors(vph))

"""
	vph_to_rs!(vph)

Convert `vph` to range and slow-time data.
"""
function vph_to_rs!(vph::VPH)
    if !vph.range_domain
        vph.data[:, :] = vph.kernel_sign < 0 ? sifft(vph.data, 1) : sfft(vph.data, 1)
        vph.range_domain = true
    end
    return nothing
end

"""
	vph_to_rs(vph)
"""
function vph_to_rs(vph::VPH)
    rsdata = deepcopy(vph)
    vph_to_rs!(rsdata)
    return rsdata
end

"""
	rs_to_vph!(vph)

Convert range and slow-time data back to frequency and slow-time data.
"""
function rs_to_vph!(vph::VPH)
    if vph.range_domain
        vph.data[:, :] = vph.kernel_sign < 0 ? sfft(vph.data, 1) : sifft(vph.data, 1)
        vph.range_domain = false
    end
    return nothing
end

"""
	rs_to_vph(vph)

Convert range and slow-time data back to frequency and slow-time data.
"""
function rs_to_vph(vph::VPH)
    new_vph = deepcopy(vph)
    rs_to_vph!(new_vph)
    return new_vph
end

"""
	freq_truncate!(vph,keep_idxs)

Truncate frequency support of VPH to range of indices.
"""
function freq_truncate!(vph::VPH, keep_idxs::UnitRange{Int})
    # confirm we are in frequency domain
    rs_to_vph!(vph)
    vph.data = vph.data[keep_idxs, :]
    vph.freq_list_hz = vph.freq_list_hz[keep_idxs]
    return nothing
end

"""
    freq_truncate!(vph, start_freq_hz, stop_freq_hz)

Truncate frequency support of VPH to the range of frequencies specified by `start_freq_hz` and `stop_freq_hz`.
"""
function freq_truncate!(vph::VPH, start_freq_hz::Real, stop_freq_hz::Real)
    start_idx = findmin(x -> abs(x - start_freq_hz), vph.freq_list_hz)[2]
    stop_idx = findmin(x -> abs(x - stop_freq_hz), vph.freq_list_hz)[2]
    freq_truncate!(vph, start_idx:stop_idx)
    return nothing
end

"""
	time_truncate!(vph, keep_idxs)

Truncate time support of VPH to range of pulses.
"""
function time_truncate!(vph::VPH, keep_idxs::AbstractVector{Int})
    vph.data = vph.data[:, keep_idxs]
    vph.ref_range_m = vph.ref_range_m[keep_idxs]
    vph.slow_time_s = vph.slow_time_s[keep_idxs]
    vph.tx_pos_enu = vph.tx_pos_enu[:, keep_idxs]
    vph.rx_pos_enu = vph.rx_pos_enu[:, keep_idxs]
    if size(vph.srp_enu, 2) > 1
        vph.srp_enu = vph.srp_enu[:, keep_idxs]
    end
    return nothing
end

"""
    time_truncate!(vph, start_time_s, stop_time_s)

Truncate VPH in time, as indicated by start and stop times.
"""
function time_truncate!(vph::VPH, start_time_s::Real, stop_time_s::Real)
    start_idx = findmin(x -> abs(x - start_time_s), vph.slow_time_s)[2]
    stop_idx = findmin(x -> abs(x - stop_time_s), vph.slow_time_s)[2]
    time_truncate!(vph, start_idx:stop_idx)
    return nothing
end

"""
    freq_pad!(vph, num_top, num_bottom = num_top)

Zero-pad VPH in frequency domain, adding `num_top` frequency bins at the low end
and `num_bottom` frequency bins at the high end.
"""
function freq_pad!(vph::VPH, num_top::Int, num_bottom::Int = num_top)
    # Ensure we are working in frequency domain
    rs_to_vph!(vph)

    # Compute frequency step
    freq_step = vph.freq_list_hz[2] - vph.freq_list_hz[1]

    # Pad data with zeros
    vph.data = vcat(
        zeros(eltype(vph.data), num_top, size(vph.data, 2)),
        vph.data,
        zeros(eltype(vph.data), num_bottom, size(vph.data, 2))
    )

    # Extend frequency list
    neg_freqs = vph.freq_list_hz[1] .+ ((-num_top):-1) .* freq_step
    pos_freqs = vph.freq_list_hz[end] .+ (1:num_bottom) .* freq_step
    vph.freq_list_hz = vcat(neg_freqs, vph.freq_list_hz, pos_freqs)
    return nothing
end

"""
    range_pad!(vph, num_top, num_bottom = num_top)

Zero-pad VPH in range (fast-time) domain, adding `num_top` range bins at the top
and `num_bottom` range bins at the bottom.
"""
function range_pad!(vph::VPH, num_top::Int, num_bottom::Int = num_top)
    # Ensure we are working in range domain
    vph_to_rs!(vph)

    # Pad data with zeros in range dimension
    vph.data = vcat(
        zeros(eltype(vph.data), num_top, size(vph.data, 2)),
        vph.data,
        zeros(eltype(vph.data), num_bottom, size(vph.data, 2))
    )

    # Update frequency list
    num_freqs_pad = size(vph.data, 1)
    df_pad = bandwidth(vph) / num_freqs_pad
    vph.freq_list_hz = center_freq(vph) .+ fft_spacing(num_freqs_pad) .* df_pad
    return nothing
end

"""
    extract_cpi(vph, cpi_idxs)

Extract a sub-VPH from the given pulse indices.
This is sometimes referred to as a Coherent Processing Interval (CPI).
"""
function extract_cpi(vph::VPH, cpi_idxs::AbstractVector{Int})
    srp_enu = size(vph.srp_enu, 2) > 1 ? vph.srp_enu[:, cpi_idxs] : vph.srp_enu
    return VPH(
        vph.data[:, cpi_idxs],
        vph.freq_list_hz,
        vph.ref_range_m[cpi_idxs],
        vph.slow_time_s[cpi_idxs],
        vph.tx_pos_enu[:, cpi_idxs];
        ref_time_s = vph.ref_time_s,
        rx_pos_enu = vph.rx_pos_enu[:, cpi_idxs],
        srp_enu = srp_enu,
        srp_lla = vph.srp_lla,
        kernel_sign = vph.kernel_sign,
        c_eff_ms = vph.c_eff_ms,
        range_domain = vph.range_domain
    )
end

"""
	show(io, vph)

Custom display method for VPH.
"""
function Base.show(io::IO, vph::VPH{T}) where {T <: AbstractFloat}
    println(io, size(vph, 1), "x", size(vph, 2), " VPH{", T, "}:")
    println(
        io,
        "    data: ",
        size(vph.data),
        @sprintf(" [%.2f%+.2fim %.2f%+.2fim ... %.2f%+.2fim %.2f%+.2fim]",
            real(vph.data[1]),
            imag(vph.data[1]),
            real(vph.data[2]),
            imag(vph.data[2]),
            real(vph.data[end - 1]),
            imag(vph.data[end - 1]),
            real(vph.data[end]),
            imag(vph.data[end]))
    )
    println(
        io,
        "    freq_list_hz: ",
        size(vph.freq_list_hz),
        @sprintf(" [%.2f ... %.2f]", vph.freq_list_hz[1], vph.freq_list_hz[end])
    )
    println(
        io,
        "    ref_range_m: ",
        size(vph.ref_range_m),
        @sprintf(" [%.2f ... %.2f]", vph.ref_range_m[1], vph.ref_range_m[end])
    )
    println(
        io,
        "    slow_time_s: ",
        size(vph.slow_time_s),
        @sprintf(" [%.2f ... %.2f]", vph.slow_time_s[1], vph.slow_time_s[end])
    )
    println(io, "    ref_time_s: ", vph.ref_time_s)
    println(
        io,
        "    tx_pos_enu: ",
        size(vph.tx_pos_enu),
        @sprintf(" [%.2f %.2f %.2f] ... [%.2f %.2f %.2f]",
            vph.tx_pos_enu[1, 1],
            vph.tx_pos_enu[2, 1],
            vph.tx_pos_enu[3, 1],
            vph.tx_pos_enu[1, end],
            vph.tx_pos_enu[2, end],
            vph.tx_pos_enu[3, end])
    )
    println(
        io,
        "    rx_pos_enu: ",
        size(vph.rx_pos_enu),
        @sprintf(" [%.2f %.2f %.2f] ... [%.2f %.2f %.2f]",
            vph.rx_pos_enu[1, 1],
            vph.rx_pos_enu[2, 1],
            vph.rx_pos_enu[3, 1],
            vph.rx_pos_enu[1, end],
            vph.rx_pos_enu[2, end],
            vph.rx_pos_enu[3, end])
    )
    if size(vph.srp_enu, 2) > 1
        println(
            io,
            "    srp_enu: ",
            size(vph.srp_enu),
            @sprintf(" [%.2f %.2f %.2f] ... [%.2f %.2f %.2f]",
                vph.srp_enu[1, 1],
                vph.srp_enu[2, 1],
                vph.srp_enu[3, 1],
                vph.srp_enu[1, end],
                vph.srp_enu[2, end],
                vph.srp_enu[3, end])
        )
    else
        println(
            io,
            "    srp_enu: ",
            size(vph.srp_enu),
            @sprintf(" [%.2f %.2f %.2f]",
                vph.srp_enu[1, 1],
                vph.srp_enu[2, 1],
                vph.srp_enu[3, 1])
        )
    end
    println(
        io,
        "    srp_lla: ",
        size(vph.srp_lla),
        @sprintf(" [%.2f %.2f %.2f]", vph.srp_lla[1], vph.srp_lla[2], vph.srp_lla[3])
    )
    println(io, "    kernel_sign: ", vph.kernel_sign)
    println(io, "    c_eff_ms: ", vph.c_eff_ms)
    println(io, "    range_domain: ", vph.range_domain)
end
