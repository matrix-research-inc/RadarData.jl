#=
	Functions to assist with forming slant plane imagery, change coordinate systems.

    License ID: SEAL_B
=#
export cart_to_horizon,
       cart_to_horizon_unwrap,
       horizon_to_cart,
       azimuth_rad,
       elevation_rad,
       slant_plane_rot_2D,
       slant_plane_rot_3D

col_norm(mat::AbstractMatrix) = sqrt.(sum(abs2, mat, dims = 1))

"""
    unit_vec(az_rad, el_rad)

Reconstruct unit vector from azimuth and elevation angle in radians.
"""
function unit_vec(az_rad, el_rad)
    return SVector(cos(el_rad) * cos(az_rad), cos(el_rad) * sin(az_rad), sin(el_rad))
end

"""
	cart_to_horizon(cart_pos)

Convert Cartesian coordinates stored as a 3 x N array to horizontal
coordinates, (range, azimuth, elevation). Azimuth is measured as
the angle with the x-axis, elevation is measured from the xy-plane.
Angles are in radians.
"""
function cart_to_horizon(cart_pos::Matrix)
    horz_pos = zeros(size(cart_pos))
    horz_pos[1, :] = col_norm(cart_pos)
    horz_pos[2, :] = atan.(cart_pos[2, :], cart_pos[1, :])
    horz_pos[3, :] = pi / 2 .- acos.(cart_pos[3, :] ./ horz_pos[1, :])
    return horz_pos
end

"""
	cart_to_horizon_unwrap(cart_pos)

Convert Cartesian coordinates stored as a 3 x N array to horizontal
coordinates, (range, azimuth, elevation). Azimuth is measured as
the angle with the x-axis, elevation is measured from the xy-plane.
Angles are in radians.
Azimuth and elevation angles are unwrapped using `DSP.unwrap`.
"""
function cart_to_horizon_unwrap(cart_pos::Matrix)
    horz_pos = cart_to_horizon(cart_pos)
    horz_pos[2, :] = unwrap(horz_pos[2, :])
    horz_pos[3, :] = unwrap(horz_pos[3, :])
    return horz_pos
end

"""
	horizon_to_cart(horz_pos)

Convert horizontal coordinates stored as a 3 x N array to Cartesian
coordinates. Angles are assumed to be in radians.
"""
function horizon_to_cart(horz_pos::Matrix)
    cart_pos = zeros(size(horz_pos))
    cart_pos[1, :] = horz_pos[1, :] .* cos.(horz_pos[2, :]) .* cos.(horz_pos[3, :])
    cart_pos[2, :] = horz_pos[1, :] .* sin.(horz_pos[2, :]) .* cos.(horz_pos[3, :])
    cart_pos[3, :] = horz_pos[1, :] .* sin.(horz_pos[3, :])
    return cart_pos
end

"""
	azimuth_rad(bisectors)

The unwrapped azimuth angle associated with the bistatic bisector stored in each column of input, in radians.
"""
azimuth_rad(bsts::Matrix) = unwrap(cart_to_horizon(bsts)[2, :])

"""
	azimuth_deg(bisectors)

The unwrapped azimuth angle associated with the bistatic bisector stored in each column of input, in degrees.
"""
azimuth_deg(bsts::Matrix) = 180.0 / pi * azimuth_rad(bsts)

"""
	elevation_rad(bisectors)

The unwrapped elevation angle associated with the bistatic bisector stored in each column of input, in radians.
"""
elevation_rad(bsts::Matrix) = unwrap(cart_to_horizon(bsts)[3, :])

"""
	elevation_deg(bisectors)

The unwrapped elevation angle associated with the bistatic bisector stored in each column of input, in degrees.
"""
elevation_deg(bsts::Matrix) = 180.0 / pi * elevation_rad(bsts)

"""
	slant_plane_rot_2D(bisectors)

Compute rotation matrix to align bisectors to form a slant-plane-aligned ground-plane image. This centers the bisectors in azimuth.
"""
function slant_plane_rot_2D(bisectors::Matrix)
    horizon_coords = cart_to_horizon(bisectors)
    az_list = unwrap(vec(horizon_coords[2, :]))
    mid_az = az_list[div(length(az_list), 2) + 1]
    return SMatrix{3, 3}(cos(mid_az), -sin(mid_az), 0, sin(mid_az), cos(mid_az), 0, 0, 0, 1)
end

"""
	slant_plane_rot_3D(bisectors, align_roll = false)

Compute rotation matrix to align bisectors to form a slant-plane aligned image.
By default, the rotation will only rotate in the azimuth and elevation directions.
To also correct for a residual roll angle between the ground plane and the slant plane, set `align_roll = true`.
"""
function slant_plane_rot_3D(bisectors::Matrix, align_roll::Bool = false)
    # Remove central azimuth and elevation components
    horizon_coords = cart_to_horizon(bisectors)
    az_list = unwrap(vec(horizon_coords[2, :]))
    el_list = unwrap(vec(horizon_coords[3, :]))
    mid_pulse = div(length(az_list), 2) + 1
    mid_az = az_list[mid_pulse]
    mid_el = el_list[mid_pulse]
    az_rot = SMatrix{3, 3}(
        cos(mid_az), -sin(mid_az), 0, sin(mid_az), cos(mid_az), 0, 0, 0, 1)
    el_rot = SMatrix{3, 3}(
        cos(mid_el), 0, -sin(mid_el), 0, 1, 0, sin(mid_el), 0, cos(mid_el))
    rot = el_rot * az_rot

    if align_roll
        # Now rotate around the x-axis to align with ground plane
        rot_bisectors = rot * bisectors
        roll_angle = atan(sum(rot_bisectors[2, :] .* rot_bisectors[3, :]) /
                          sum(abs2, rot_bisectors[2, :]))
        # refCoord = rot_bisectors[:, mid_pulse]
        # lsv, sings, rsv = svd(rot_bisectors .- refCoord)
        # roll_angle = atan(lsv[3] / lsv[2])
        roll_rot = SMatrix{3, 3}(
            1,
            0,
            0,
            0,
            cos(roll_angle),
            -sin(roll_angle),
            0,
            sin(roll_angle),
            cos(roll_angle)
        )
        rot = roll_rot * rot
    end
    return rot
end
