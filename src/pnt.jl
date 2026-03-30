#=
    Structures for managing external Position, Navigation and Timing information

    License ID: SEAL_B
=#

export PVA,
       PVAHeader,
       PNT,
       PNTHistory,
       receive_pva!,
       interp_states,
       interp_ecef,
       range_to_point,
       range_rate_to_point,
       ned_to_enu,
       body_to_ecef,
       enu_rotation,
       intersect_los_plane

"""
Position, Velocity, Attitude message header.
Corresponds to ASPN convention in `aspn_schema/header.schema.json`
"""
struct PVAHeader
    vendor_id::UInt32
    device_id::UInt64
    context_id::UInt32
    sequence_id::UInt16
end

"""
Position, Velocity, Attitude message.
Corresponds to ASPN convention in `aspn_schema/measurement_position_velocity_attitude.schema.json`.
"""
struct PVA
    header::PVAHeader
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
    attitude::SVector{3, Float64}
    covariance::SMatrix{9, 9, Float64}
    time_of_validity::Float64
end

"""
    Position, navigation, and timing data. Float64 was chosen to reduce precision errors in position and timing.

ecef = Earth-Centered, Earth-Fixed coordinate.
vel_enu_ms = Velocity of platform in local East-North-Up coordinates.
attitude = rotation of body frame relative to local NED coordinate frame.
time_s = Time of measurement, in seconds. Typically maintained relative to reference time.
"""
@kwdef struct PNT
    ecef::ECEF{Float64} = ECEF{Float64}(NaN, NaN, NaN)
    vel_enu_ms::ENU{Float64} = ENU(NaN, NaN, NaN)
    attitude::RotXYZ{Float64} = RotXYZ(0.0, 0.0, 0.0)
    time_s::Float64 = NaN
    lla::LLA{Float64} = LLA(ecef, wgs84)
    vel_ecef_ms::ECEF{Float64} = ECEF(vel_enu_ms, ecef, wgs84) - ecef
end

function PNT(pva::PVA, ref_time_s::Float64 = 0.0)
    lla = LLA(pva.position[1], pva.position[2], pva.position[3])
    return PNT(;
        ecef = ECEF(lla, wgs84),
        vel_enu_ms = ENU(pva.velocity[2], pva.velocity[1], -pva.velocity[3]),
        attitude = RotXYZ(pva.attitude[1], pva.attitude[2], pva.attitude[3]),
        time_s = pva.time_of_validity - ref_time_s,
        lla
    )
end

"""
    PNT history maintains a buffer of previous `max_states` PNT states.

`current_end` is an index indicating the last-stored state.
When `current_end` exceeds `max_states`, `current_end` is reset to 1 and the buffer begins overwriting the oldest elements.
Also includes `ref_time_s`; all contained states are assumed to have `time_s` relative to this time.
"""
@kwdef mutable struct PNTHistory
    max_states::Int = 1000
    current_end::Int = 0
    ref_time_s::Float64 = time()
    states::Vector{PNT} = fill(PNT(), max_states)
end

Base.getindex(history::PNTHistory, idx) = history.states[idx]
function Base.push!(history::PNTHistory, pnt::PNT)
    history.current_end = mod1(history.current_end + 1, history.max_states)
    history.states[history.current_end] = pnt
    return nothing
end

"""
    receive_pva!(history, pva)

Receive a position, velocity, acceleration message. Overwrites the oldest state if PNTHistory has reached its maximum size.
"""
function receive_pva!(history::PNTHistory, pva::PVA)
    push!(history, PNT(pva, history.ref_time_s))
end

"""
    get_closest_states(history, time_s)

Identify the states which are closest to the input relative time.
"""
function get_closest_states(history::PNTHistory, time_s::Float64)
    idx = history.current_end
    ref = mod1(idx + 1, history.max_states)
    # If the requested time is too recent
    if history[idx].time_s < time_s
        return history[idx], nothing
    end
    # Otherwise back-track until we get a time before the requested time
    while (history[idx].time_s > time_s)
        idx = mod1(idx - 1, history.max_states)
        if idx == ref
            # Prevent never-ending loop
            return nothing, nothing
        end
    end
    return history[idx], history[mod1(idx + 1, history.max_states)]
end

"""
    get_states(history, time_start, time_end)

Identify the indices of the states which cover the input time range.
"""
function get_states(history::PNTHistory, time_start::Float64, time_end::Float64)
    idx_e = history.current_end
    ref = mod1(idx_e + 1, history.max_states)
    # Back-track until we get a time before the requested time
    while history[idx_e].time_s > time_end
        idx_e = mod1(idx_e - 1, history.max_states)
        if idx_e == ref
            # Prevent never-ending loop, end time predates buffer
            return nothing, nothing
        end
    end
    idx_s = copy(idx_e)
    if idx_e != history.current_end
        # Need to increment by 1 to include the end time
        idx_e = mod1(idx_e + 1, history.max_states)
    end
    # Keep going to capture start time
    while history[idx_s].time_s > time_start
        idx_s = mod1(idx_s - 1, history.max_states)
        if idx_s == ref
            # Prevent never-ending loop, start time predates buffer
            break
        end
    end
    return idx_s, idx_e
end

"""
    interp_states(history, time_s)

Return best linearly interpolated estimate of state at the relative input time, in seconds.
Assumes `time_s` is a single float or is sorted in increasing order.
No extrapolation.
"""
function interp_states(history::PNTHistory, time_s::Float64)
    left_state, right_state = get_closest_states(history, time_s)
    if isnothing(left_state)
        # Time was too old
        return nothing
    elseif isnothing(right_state)
        # Time was too new, don't extrapolate
        return left_state
    else
        # Time was in between two states, linearly interp
        interval_s = right_state.time_s - left_state.time_s
        rc = (time_s - left_state.time_s) / interval_s
        lc = 1.0 - rc
        return PNT(;
            ecef = lc .* left_state.ecef .+ rc .* right_state.ecef,
            vel_enu_ms = lc .* left_state.vel_enu_ms .+ rc .* right_state.vel_enu_ms,
            attitude = RotXYZ(lc .* QuatRotation(left_state.attitude) .+
                              rc .* QuatRotation(right_state.attitude)),
            time_s
        )
    end
end

"""
    interp_ecef(history, time_s)

Interpolate ECEF field of position history for multiple times.
"""
function interp_ecef(history::PNTHistory, time_s::Vector{Float64})
    # Get indices which include time_s
    idx_s, idx_e = get_states(history, time_s[1], time_s[end])
    if isnothing(idx_e)
        @error "PNT history does not overlap with requested time extent."
    end

    # Collect the data, order it
    num_knots = idx_s < idx_e ? idx_e - idx_s + 1 : history.max_states - idx_s + 1 + idx_e
    knot_times = zeros(Float64, num_knots)
    knot_vals = zeros(Float64, 3, num_knots)
    if idx_s < idx_e
        knot_times[:] = map(x -> x.time_s, history.states[idx_s:idx_e])
        knot_vals[:, :] = stack(x -> x.ecef, history.states[idx_s:idx_e])
    else
        knot_times[1:(history.max_states - idx_s + 1)] = map(
            x -> x.time_s, history.states[idx_s:end])
        knot_times[(history.max_states - idx_s + 2):end] = map(
            x -> x.time_s, history.states[1:idx_e])
        knot_vals[:, 1:(history.max_states - idx_s + 1)] = stack(
            x -> x.ecef, history.states[idx_s:end])
        knot_vals[:, (history.max_states - idx_s + 2):end] = stack(
            x -> x.ecef, history.states[1:idx_e])
    end

    # Now interpolate
    itp = extrapolate(
        interpolate((1:3, knot_times), knot_vals, (NoInterp(), Gridded(Linear()))),
        NaN
    )
    out_ecef = zeros(3, length(time_s))
    for tidx in eachindex(time_s)
        out_ecef[1, tidx] = itp(1, time_s[tidx])
        out_ecef[2, tidx] = itp(2, time_s[tidx])
        out_ecef[3, tidx] = itp(3, time_s[tidx])
    end
    return out_ecef
end

"""
    range_to_point(pnt,ecef)

Compute one-way range (m) from PNT to ECEF coordinate.
"""
function range_to_point(pnt::PNT, ecef::ECEF{Float64})
    return sqrt((pnt.ecef[1] - ecef[1])^2 +
                (pnt.ecef[2] - ecef[2])^2 +
                (pnt.ecef[3] - ecef[3])^2)
end

"""
    range_rate_to_point(pnt, ecef)

Compute one-way range (m) and one-way change in range (m/s) from PNT to ECEF coordinate.
"""
function range_rate_to_point(pnt::PNT, ecef::ECEF{Float64})
    vel_ecef_ms = pnt.vel_ecef_ms
    range_to_pt = range_to_point(pnt, ecef)
    vel_to_pt = ((pnt.ecef[1] - ecef[1]) * vel_ecef_ms[1] +
                 (pnt.ecef[2] - ecef[2]) * vel_ecef_ms[2] +
                 (pnt.ecef[3] - ecef[3]) * vel_ecef_ms[3]) / range_to_pt
    return range_to_pt, vel_to_pt
end

"""
    range_rate_to_point(pnt, enu, ecef_ref)

Compute one-way range (m) and one-way change in range (m/s) from PNT to ENU coordinate.
"""
function range_rate_to_point(
        pnt::PNT,
        enu::ENU{Float64},
        ecef_ref::ECEF{Float64} = pnt.ecef
)
    enu_pos = ENU(pnt.ecef, ecef_ref, wgs84)
    range_to_pt = sqrt((enu_pos[1] - enu[1])^2 + (enu_pos[2] - enu[2])^2 +
                       (enu_pos[3] - enu[3])^2)
    vel_to_pt = ((enu_pos[1] - enu[1]) * pnt.vel_enu_ms[1] +
                 (enu_pos[2] - enu[2]) * pnt.vel_enu_ms[2] +
                 (enu_pos[3] - enu[3]) * pnt.vel_enu_ms[3]) / range_to_pt
    return range_to_pt, vel_to_pt
end

"""
    ned_to_enu(ned)

Convert North-East-Down coordinates to East-North-Up.
"""
function ned_to_enu(ned::V) where {V <: AbstractVector{Float64}}
    return ENU(ned[2], ned[1], -ned[3])
end

"""
    body_to_ecef(ned)

Convert body coordinates to ECEF orientations.
"""
function body_to_ecef(body::V, pnt::PNT) where {V <: AbstractVector{Float64}}
    ned = pnt.attitude * body
    return ECEF(ned_to_enu(ned), pnt.ecef, wgs84) - pnt.ecef
end

"""
    enu_rotation(lla)

Compute rotation matrix from local ENU coordinates to ECEF coordinates, given `lla` as origin of ENU frame.
`ecef = R * enu + ecef_ref`
ENU coordinate axes are described by the columns of `R`.
"""
function enu_rotation(lla::LLA{Float64})
    return RotZX(-(lla.lon + 90) * pi / 180, (lla.lat - 90) * pi / 180)
end

function enu_rotation(ecef::ECEF{Float64})
    return enu_rotation(LLA(ecef, wgs84))
end

"""
    intersect_ray_wgs84(origin, direction)

Compute intersection of ray in ECEF with WGS84 ellipsoid.
"""
function intersect_ray_wgs84(origin::ECEF, direction::ECEF)
    o_x, o_y, o_z = origin
    d_x, d_y, d_z = direction
    a_e = wgs84_ellipsoid.a
    b_e = wgs84_ellipsoid.b

    # Intersection of line with ellipsoid results in quadratic equation to solve
    aquad = (d_x^2 / a_e^2) + (d_y^2 / b_e^2) + (d_z^2 / b_e^2)
    bquad = 2 * ((o_x * d_x / a_e^2) + (o_y * d_y / b_e^2) + (o_z * d_z / b_e^2))
    cquad = ((o_x^2 / a_e^2) + (o_y^2 / b_e^2) + (o_z^2 / b_e^2)) - 1
    discriminant = bquad^2 - 4 * aquad * cquad

    if discriminant < 0
        return NaN
    else
        t1 = (-bquad - sqrt(discriminant)) / (2 * aquad)
        t2 = (-bquad + sqrt(discriminant)) / (2 * aquad)
        return origin + direction * min(abs(t1), abs(t2))
    end
end

"""
    intersect_los_plane(body_los, pnt, ground_hae_m = 0)

Compute intersection of line of sight vector relative to body frame with ENU plane
referenced to `pnt.ecef`, and return ENU coordinate. Ground plane is assumed to have
given height above ellipsoid (HAE).
"""
function intersect_los_plane(
        body_los::V,
        pnt::PNT,
        ground_hae_m::Float64 = 0.0
) where {V <: AbstractVector{Float64}}
    # convert body vector to ENU frame
    los_enu = ned_to_enu(pnt.attitude * body_los)
    # compute intersection with ground plane, assuming given HAE
    if los_enu[3] < 0
        ray_len = (pnt.lla.alt - ground_hae_m) / abs(los_enu[3])
        return ray_len * los_enu
    else
        return NaN * los_enu
    end
end

# """
#     PNTHistory(file_name)

# Load PNT data from file.
# """
# function PNTHistory(file_name::String)
#     _, ext = splitext(file_name)
#     if ext == ".h5"
#         fid = h5open(file_name, "r")
#         attitude_rad = read(fid["attitude_rad"])
#         ecf_coords_m = read(fid["ecf_coords_m"])
#         vel_ned_ms = read(fid["vel_ned_ms"])
#         time_s = read(fid["time_s"])
#         close(fid)

#         history = PNTHistory(length(time_s))
#         history.current_end = history.max_states
#         for idx in 1:history.max_states
#             history.states[idx] = PNT(;
#                 ecef = ECEF{Float64}(ecf_coords_m[:, idx]),
#                 vel_enu_ms = SVector{3,Float64}(
#                     vel_ned_ms[2, idx],
#                     vel_ned_ms[1, idx],
#                     -vel_ned_ms[3, idx]
#                 ),
#                 attitude = RotXYZ{Float64}(
#                     attitude_rad[1, idx],
#                     attitude_rad[2, idx],
#                     attitude_rad[3, idx]
#                 ),
#                 time_s = Float64(time_s[idx])
#             )
#         end
#     else
#         error("File format not supported.")
#     end
#     return history
# end
