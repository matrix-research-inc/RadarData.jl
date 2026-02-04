export push_status!, push_waypoint!

using Geodesy: ENU, ECEF, wgs84
using Rotations: RotXYZ

function gps_to_unixtime(gps_time_s::Float64)
    gps_epoch_s = 3.159648e8
    # Add 18 leap seconds, should work for a while
    return gps_epoch_s + gps_time_s + 18
end

"""
    push_status!(history, status)
    
Push status data from an `INSFeed` to a `PNTHistory`.
"""
function push_status!(history::PNTHistory, status)
    ecef = ECEF{Float64}(status.gnss.gnss1PosEcef)
    vel_enu_ms = ENU{Float64}(
        status.gnss.gnss1VelNed[2],
        status.gnss.gnss1VelNed[1],
        -status.gnss.gnss1VelNed[3]
    )
    utc_time = status.gnss.gnss1TimeUtc
    # Year is signed output from the year 2000
    date_time = DateTime(
        2000 + Int64(utc_time.year),
        utc_time.month,
        utc_time.day,
        utc_time.hour,
        utc_time.minute,
        utc_time.second
    )
    # fracSec are milliseconds
    time_s =
        datetime2unix(date_time) + utc_time.fracSec / 1000 - history.ref_time_s
    push!(history, PNT(; ecef, vel_enu_ms, time_s))
    return nothing
end

"""
    push_waypoint!(history, waypoint)
    
Push waypoint data from an `INSFeed` to a `PNTHistory`.
"""
function push_waypoint!(history::PNTHistory, waypoint)
    ecef = ECEF{Float64}(waypoint.ins.posEcef)
    vel_ecef_ms = ECEF{Float64}(waypoint.ins.velEcef)
    vel_enu_ms = ENU(vel_ecef_ms .+ ecef, ecef, wgs84)
    # GPS time comes as nanoseconds
    gps_time_s = Float64(waypoint.time.timeGps) * 1e-9
    time_s = gps_to_unixtime(gps_time_s) - history.ref_time_s
    attitude = RotXYZ(
        waypoint.attitude.ypr[3] * pi / 180,
        waypoint.attitude.ypr[2] * pi / 180,
        waypoint.attitude.ypr[1] * pi / 180
    )
    push!(history, PNT(; ecef, vel_enu_ms, time_s, attitude))
    return nothing
end