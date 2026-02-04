"""An extension to RadarData for loading data saved for the DARPA Moving Target Recognition (MTR) program into VPH format."""
module MTR

using MAT
using Geodesy
using RadarData: VPH

"""
	read_vph_legacy(file_name)

Read older version of video phase history array from file. Used for compatibility with MTR datasets and MATLAB repository.
"""
function read_vph_legacy(file_name::String)
    fid = open(file_name, "r")

    # Get length of vph (number of channels)
    num_vph = read(fid, Int8)
    vph_array = Array{VPH,1}(undef, num_vph)

    for ii in 1:N
        # Get size of data
        n_freqs = read(fid, Int64)
        n_pulses = read(fid, Int64)

        vph = VPH(n_freqs, n_pulses)
        read!(fid, vph.data)
        read!(fid, vph.freq_list_hz)
        read!(fid, vph.ref_range_m)
        vph.c_eff_ms = read(fid, Float64)
        read!(fid, vph.slow_time_s)
        vph.ref_time_s = read(fid, Float64)
        read!(fid, vph.tx_pos_enu)
        read!(fid, vph.rx_pos_enu)
        vph.kernel_sign = read(fid, Int8)

        vph_array[ii] = vph

        # Remaining fields for backwards compat.
        pri = read(fid, Float64)
        pulse_idxs = Array{Float64,1}(undef, n_pulses)
        read!(fid, pulse_idxs)
        slant_plane_rot_2D = Array{Float64,2}(undef, 3, 3)
        read!(fid, slant_plane_rot_2D)
        A = read(fid, Int8)
        slant_plane_rot_3D = Array{Float64,2}(undef, 3, 3)
        read!(fid, slant_plane_rot_3D)
    end
    close(fid)
    return vph_array
end

"""
    load_mtr_unpfa(file_name)

Load "unpfa-ed" data from the DARPA Moving Target Recognition (MTR) program.
Return data is in RadarData.VPH format.
"""
function load_mtr_unpfa(file_name::String)
    vph_mat = matread(file_name)
    freq_list_hz =
        vph_mat["ph_meta_unPFA"]["frequency_start"][1] .+
        (0:size(vph_mat["ph_data_unPFA"]["data"], 1)-1) .*
        vph_mat["ph_meta_unPFA"]["frequency_step"][1]
    slow_time_s =
        (0:size(vph_mat["ph_data_unPFA"]["data"], 2)-1) .*
        vph_mat["ph_meta_unPFA"]["pulse_sample_time"][1]
    ref_range_m = 2 .* vec(vph_mat["ph_meta_unPFA"]["mocomp_range_offset"])

    # Convert ECEF to ENU
    tx_ecef = vph_mat["ph_meta_unPFA"]["tx_position"]
    rx_ecef = vph_mat["ph_meta_unPFA"]["rx_position"]
    target_ecef = vph_mat["target_truth"]["position_ecf"]
    tx_enu = similar(tx_ecef)
    rx_enu = similar(rx_ecef)
    for idx in axes(tx_ecef, 2)
        tx_enu[:, idx] =
            ENU(ECEF(tx_ecef[:, idx]), ECEF(target_ecef[:, idx]), wgs84)
        rx_enu[:, idx] =
            ENU(ECEF(rx_ecef[:, idx]), ECEF(target_ecef[:, idx]), wgs84)
    end

    return VPH(
        vph_mat["ph_data_unPFA"]["data"],
        collect(freq_list_hz),
        ref_range_m,
        slow_time_s,
        tx_enu;
        ref_time_s = vph_mat["ph_meta_unPFA"]["pulse_start_time"][1],
        rx_pos_enu = rx_enu,
        kernel_sign = -1.0,
        c_eff_ms = vph_mat["ph_meta_unPFA"]["c_eff"]
    )
end

end