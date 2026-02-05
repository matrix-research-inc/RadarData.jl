#=
	Methods for motion compensating video phase history

    License ID: SEAL_B
=#
export phase_to_range, range_to_phase, motion_compensate!, motion_compensate

"""
    phase_to_range(phase, center_freq, c_eff_ms)

Convert phase in radians to range in meters.
"""
phase_to_range(phase::Vector, center_freq::Real, c_eff_ms::Real) = (c_eff_ms / center_freq /
                                                                    2 / pi) .* phase

"""
    phase_to_range(phase, vph)

Convert phase in radians to range in meters.
"""
phase_to_range(phase::Vector, vph::VPH) = phase_to_range(
    phase, center_freq(vph), vph.c_eff_ms)

"""
    range_to_phase(phase, center_freq, c_eff_ms)

Convert range in meters to phase in radians.
"""
range_to_phase(range::Vector, center_freq::Real, c_eff_ms::Real) = (2 * pi * center_freq /
                                                                    c_eff_ms) .* range

"""
    range_to_phase(phase, vph)

Convert range in meters to phase in radians.
"""
range_to_phase(range::Vector, vph::VPH) = range_to_phase(
    range, center_freq(vph), vph.c_eff_ms)

"""
	motion_compensate(data_matrix, k_list, range_list; kernel_sign = -1)

Apply circular shift motion compensation directly to phase history data matrix.
"""
function motion_compensate!(
        data::Matrix,
        k_list::Vector{Float64},
        range_list::Vector{Float64};
        kernel_sign = -1
)
    data[:, :] .*= cis.(kernel_sign .* k_list .* range_list')
    return nothing
end

"""
	motion_compensate!(vph, k_list, range_list)

Multiply `vph` by `exp( -ikr )` where `k` is drawn from `k_list` and `r` is drawn from `range_list`.
This implements a circular shift of `-r` in the range direction of the range / slowtime data corresponding to `vph`.
Note that `range_list` is the two-way range, so users correcting monostatic data will need to account for this in `range_list`.
"""
function motion_compensate!(vph::VPH, k_list::Vector{Float64}, range_list::Vector{Float64})
    was_rs = vph.range_domain
    rs_to_vph!(vph)
    motion_compensate!(vph.data, k_list, range_list; kernel_sign = vph.kernel_sign)
    vph.ref_range_m[:] .-= range_list
    if was_rs
        vph_to_rs!(vph)
    end
    return nothing
end

"""
	motion_compensate(vph, k_list, range_list)
"""
function motion_compensate(vph::VPH, k_list::Vector{Float64}, range_list::Vector{Float64})
    new_vph = deepcopy(vph)
    motion_compensate!(new_vph, k_list, range_list)
    return new_vph
end

"""
	motion_compensate!(vph, range_list)
"""
function motion_compensate!(vph::VPH, range_list::Vector{Float64})
    k_list = spatial_freqs(vph)
    motion_compensate!(vph, k_list, range_list)
    return nothing
end

"""
	motion_compensate(vph, range_list)
"""
function motion_compensate(vph::VPH, range_list::Vector{Float64})
    k_list = spatial_freqs(vph)
    return motion_compensate(vph, k_list, range_list)
end
