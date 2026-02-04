#=
    License ID: SEAL_B
=#

using Test
using JSON3
using Geodesy
using Rotations

# Initialize
max_states = 100
history = PNTHistory(max_states)

# Pack a test PVA message
pva_dict = Dict(
    "header" =>
        Dict("vendor_id" => 0, "device_id" => 0, "context_id" => 0, "sequence_id" => 0),
    "time_of_validity" => 0.0,
    "position" => [39.729866; -84.077144; 0],
    "velocity" => zeros(3),
    "attitude" => zeros(3),
    "covariance" => vec(zeros(9, 9))
)

JSON3.write("test_pva.json", pva_dict)
pva = JSON3.read("test_pva.json", PVA)
rm("test_pva.json")

# Push PNT state onto history 
receive_pva!(history, pva)

@test history[1].vel_enu_ms == zeros(3)

# Test interpolation of states 
history.states[1] = PNT(;
    ecef = ECEF(0, 0, 0),
    vel_enu_ms = zeros(3),
    attitude = RotXYZ(0, 0, 0),
    time_s = 0
)
history.states[2] = PNT(;
    ecef = ECEF(1, 2, 3),
    vel_enu_ms = ones(3),
    attitude = RotXYZ(0, 0, pi / 90),
    time_s = 1
)
history.current_end = 2

time_s = 0.3
itp_pnt = interp_states(history, time_s)

@test isapprox(Vector(itp_pnt.ecef), [time_s, 2 * time_s, 3 * time_s])
@test isapprox(itp_pnt.vel_enu_ms, time_s * ones(3))
@test isapprox(itp_pnt.attitude, RotXYZ(0, 0, time_s * pi / 90); rtol = 1e-6)
@test time_s == itp_pnt.time_s

# Test range computations
ecef_ref = ECEF(LLA(0, 0, 0), wgs84)
pnt = PNT(;
    ecef = ecef_ref + [1000, 0, 0],
    vel_enu_ms = SVector{3,Float64}(0, 0, -10),
    attitude = RotXYZ(0, 0, 0),
    time_s = 0
)
meas_pt_enu = ENU(0.0, 0.0, 0.0)
meas_pt_ecef = ECEF(meas_pt_enu, ecef_ref, wgs84)

range_to_pt, rate_to_pt = range_rate_to_point(pnt, meas_pt_ecef)
range_to_pt2, rate_to_pt2 = range_rate_to_point(pnt, meas_pt_enu, ecef_ref)

@test isapprox(range_to_pt, 1000)
@test isapprox(rate_to_pt, -10)
@test isapprox(range_to_pt, range_to_pt2)
@test isapprox(rate_to_pt, rate_to_pt2)

# # Load a PNT History from file 
# file_name = string(ENV["DATA_DIR"], "/alt-pnt/flight_paths/low_agl_flight.h5")
# history = PNTHistory(file_name)
