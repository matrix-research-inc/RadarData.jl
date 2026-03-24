#=
    Contains image metadata structs.
=#
export ImagePlane,
       GROUND_PLANE,
       FIXED_EL_PLANE,
       SLANT_PLANE,
       PFAKnots

"""
    ImagePlane

Enum specifying the image focus plane convention.

Values:
- `GROUND_PLANE`: Image plane coincides with ENU ground plane (z=0)
- `FIXED_EL_PLANE`: Slant plane for constant-elevation flight through aperture center
- `SLANT_PLANE`: Best-fit plane to VPH bisector vectors
"""
@enum ImagePlane GROUND_PLANE FIXED_EL_PLANE SLANT_PLANE

"""
Polar Format Algorithm Knots data structure.

This structure holds metadata relevant to the formation and inversion of polar format images.

$(TYPEDFIELDS)
"""
@kwdef mutable struct PFAKnots
    "Frequency list, radians/meter"
    spatial_freqs::Vector{Float64}
    "Bistatic bisector vectors (3 x num_pulses)"
    bisectors::Matrix{Float64}
    "Two-way reference range per pulse, in meters"
    ref_range_m::Vector{Float64}
    "Time relative to reference time, per pulse, in seconds"
    slow_time_s::Vector{Float64}
    "Reference time, in seconds since Unix Epoch"
    ref_time_s::Float64
    "Spatial frequency center, removed for basebanding (radians/meter)"
    k_center::SVector{2, Float64}
    "Frequency and slow-time indices which correspond to k-center"
    exp_inds::SVector{2, Float64}
    "Spatial frequency in first dimension, uniformly spaced for imaging (radians/meter)"
    kx_image::Vector{Float64}
    "Spatial frequency in second dimension, uniformly spaced for imaging (radians/meter)"
    ky_image::Vector{Float64}
    "Image axis first dimension, meters. Range for slant plane aligned imagery."
    x_image::Vector{Float64}
    "Image axis second dimension, meters. Cross range for slant plane aligned imagery."
    y_image::Vector{Float64}
    "Rotation from ENU to slant plane"
    slant_plane_rot::SMatrix{3, 3, Float64} = SMatrix{3, 3}(diagm(ones(3)))
    "Polynomial coefficients for slant-plane-aligned transmitter position"
    tx_poly::SMatrix{3, 3, Float64} = SMatrix{3, 3, Float64}(zeros(3, 3))
    "FFT sign convention, default is -1"
    kernel_sign::Float64 = -1.0
    "Whether to use separable interpolation or joint interpolation"
    separable::Bool = true
    "Whether to use an inner box in frequency domain or outer box"
    inner_box::Bool = false
    "Whether to apply quadratic phase error correction"
    quadratic::Bool = false
    "Whether to correct the characteristic PFA distortion by resampling the image"
    distortion::Bool = false
    "Choice of focus plane for image, see ImagePlane enum for options"
    image_plane::ImagePlane = SLANT_PLANE
    "Resamples to align image with VPH sensor coordinate convention (ENU)"
    enu_align::Bool = false
    "Derived quantity, range, rate, etc. of the sensor at aperture center"
    comp_poly::SVector{3, Float64} = compute_comp_poly(tx_poly)
    "Derived quantity, starting index for radial interpolation"
    radius_start::Vector{Float64} = (kx_image[1] .- spatial_freqs[1] .* bisectors[1, :] .+
                                     k_center[1]) ./
                                    ((spatial_freqs[2] - spatial_freqs[1]) .*
                                     bisectors[1, :]) .+ 1
    "Derived quantity, index step for radial interpolation"
    radius_delta::Vector{Float64} = (kx_image[2] - kx_image[1]) ./
                                    ((spatial_freqs[2] - spatial_freqs[1]) .*
                                     bisectors[1, :])
end

function Base.getproperty(pfa_knots::PFAKnots, name::Symbol)
    if name == :num_x
        return length(getfield(pfa_knots, :x_image))
    elseif name == :num_y
        return length(getfield(pfa_knots, :y_image))
    elseif name == :delta_kx
        return getfield(pfa_knots, :kx_image)[2] - getfield(pfa_knots, :kx_image)[1]
    elseif name == :delta_ky
        return getfield(pfa_knots, :ky_image)[2] - getfield(pfa_knots, :ky_image)[1]
    elseif name == :bandwidth_kx
        return getfield(pfa_knots, :kx_image)[end] - getfield(pfa_knots, :kx_image)[1]
    elseif name == :bandwidth_ky
        return getfield(pfa_knots, :ky_image)[end] - getfield(pfa_knots, :ky_image)[1]
    elseif name == :delta_x
        return getfield(pfa_knots, :x_image)[2] - getfield(pfa_knots, :x_image)[1]
    elseif name == :delta_y
        return getfield(pfa_knots, :y_image)[2] - getfield(pfa_knots, :y_image)[1]
    elseif name == :unit_vec_x
        return -getfield(pfa_knots, :slant_plane_rot)[1, :]
    elseif name == :unit_vec_y
        return -getfield(pfa_knots, :slant_plane_rot)[2, :]
    elseif name == :center_freq
        return getfield(pfa_knots, :spatial_freqs)[1] +
               getfield(pfa_knots, :exp_inds)[1] * (getfield(pfa_knots, :spatial_freqs)[2] -
                getfield(pfa_knots, :spatial_freqs)[1])
    else
        return getfield(pfa_knots, name)
    end
end

"""
    dot3(array, vector, colidx=1, transpose=false)

Loop-unrolled, non-heap allocating version of dot product for 3-vectors.
First argument may be an array, where `colidx` is the index of the column to use in dot product.
`transpose` will transpose the array argument prior to multiplication.
No checks in place for attempting to transpose a vector in the first argument.
"""
function dot3(
        array::AbstractArray,
        vector::AbstractVector,
        colidx::Int = 1,
        transpose::Bool = false
)
    if transpose
        return array[colidx, 1] * vector[1] +
               array[colidx, 2] * vector[2] +
               array[colidx, 3] * vector[3]
    else
        return array[1, colidx] * vector[1] +
               array[2, colidx] * vector[2] +
               array[3, colidx] * vector[3]
    end
end

"""
    compute_comp_poly(tx_poly::SMatrix{3,3,Float64}) -> SVector{3,Float64}

Compute motion compensation polynomial [r₀, ṙ₀, r̈₀] from transmitter path polynomial coefficients.

# Arguments
- `tx_poly`: 3×3 matrix of polynomial coefficients for [x(τ), y(τ), z(τ)]

# Returns
Range, range rate, and range acceleration at aperture center.
"""
function compute_comp_poly(tx_poly::SMatrix{3, 3, Float64})
    r0 = sqrt.(sum(abs2, tx_poly[:, 1]))
    dr0 = dot3(tx_poly, tx_poly[:, 2], 1) / r0
    ddr0 = (sum(abs2, tx_poly[:, 2]) + 2 * dot3(tx_poly, tx_poly[:, 3], 1) - dr0^2) / r0
    return SVector{3, Float64}(r0, dr0, ddr0)
end

"""
    Base.show(io::IO, pfa_knots::PFAKnots)

Pretty-print PFAKnots structure with formatted field display.
"""
function Base.show(io::IO, pfa_knots::PFAKnots)
    println("PFAKnots Object")
    println(
        io,
        "   VPH size: (",
        length(pfa_knots.spatial_freqs),
        ", ",
        size(pfa_knots.bisectors, 2),
        ")"
    )
    println(io, "   Image size: (", pfa_knots.num_x, ", ", pfa_knots.num_y, ")")
    println(
        io,
        "   kx_image: ",
        @sprintf(" [%.2f, %.2f ... %.2f]",
            pfa_knots.kx_image[1],
            pfa_knots.kx_image[2],
            pfa_knots.kx_image[end])
    )
    println(
        io,
        "   ky_image: ",
        @sprintf(" [%.2f, %.2f ... %.2f]",
            pfa_knots.ky_image[1],
            pfa_knots.ky_image[2],
            pfa_knots.ky_image[end])
    )
    println(
        io,
        "   x_image: ",
        @sprintf(" [%.2f, %.2f ... %.2f]",
            pfa_knots.x_image[1],
            pfa_knots.x_image[2],
            pfa_knots.x_image[end])
    )
    println(
        io,
        "   y_image: ",
        @sprintf(" [%.2f, %.2f ... %.2f]",
            pfa_knots.y_image[1],
            pfa_knots.y_image[2],
            pfa_knots.y_image[end])
    )
    println(
        io,
        "   spatial_freqs: ",
        size(pfa_knots.spatial_freqs),
        @sprintf(" [%.2f ... %.2f]",
            pfa_knots.spatial_freqs[1],
            pfa_knots.spatial_freqs[end])
    )
    println(
        io,
        "   ref_range_m: ",
        size(pfa_knots.ref_range_m),
        @sprintf(" [%.2f ... %.2f]", pfa_knots.ref_range_m[1], pfa_knots.ref_range_m[end])
    )
    println(
        io,
        "   slow_time_s: ",
        size(pfa_knots.slow_time_s),
        @sprintf(" [%.2f ... %.2f]", pfa_knots.slow_time_s[1], pfa_knots.slow_time_s[end])
    )
    println(
        io,
        "   bisectors: ",
        size(pfa_knots.bisectors),
        @sprintf(" [%.2f %.2f %.2f] ... [%.2f %.2f %.2f]",
            pfa_knots.bisectors[1, 1],
            pfa_knots.bisectors[2, 1],
            pfa_knots.bisectors[3, 1],
            pfa_knots.bisectors[1, end],
            pfa_knots.bisectors[2, end],
            pfa_knots.bisectors[3, end])
    )
    println(
        io,
        "   k_center: ",
        @sprintf(" [%.2f %.2f]", pfa_knots.k_center[1], pfa_knots.k_center[2],)
    )
    println(
        io,
        "   exp_inds: ",
        @sprintf(" [%.2f %.2f]", pfa_knots.exp_inds[1], pfa_knots.exp_inds[2],)
    )
    println(
        io,
        "   slant_plane_rot: ",
        @sprintf(" [%.2f %.2f %.2f] [%.2f %.2f %.2f] [%.2f %.2f %.2f]",
            pfa_knots.slant_plane_rot[1, 1],
            pfa_knots.slant_plane_rot[2, 1],
            pfa_knots.slant_plane_rot[3, 1],
            pfa_knots.slant_plane_rot[1, 2],
            pfa_knots.slant_plane_rot[2, 2],
            pfa_knots.slant_plane_rot[3, 2],
            pfa_knots.slant_plane_rot[1, 3],
            pfa_knots.slant_plane_rot[2, 3],
            pfa_knots.slant_plane_rot[3, 3])
    )
    println(
        io,
        "   tx_poly: ",
        @sprintf(" [%.2f %.2f %.2f] [%.2f %.2f %.2f] [%.2f %.2f %.2f]",
            pfa_knots.tx_poly[1, 1],
            pfa_knots.tx_poly[2, 1],
            pfa_knots.tx_poly[3, 1],
            pfa_knots.tx_poly[1, 2],
            pfa_knots.tx_poly[2, 2],
            pfa_knots.tx_poly[3, 2],
            pfa_knots.tx_poly[1, 3],
            pfa_knots.tx_poly[2, 3],
            pfa_knots.tx_poly[3, 3])
    )
    println(
        io,
        "   comp_poly: ",
        @sprintf("[%.2f %.2f %.2f]",
            pfa_knots.comp_poly[1],
            pfa_knots.comp_poly[2],
            pfa_knots.comp_poly[3])
    )
    println(io, "   separable: ", pfa_knots.separable)
    println(io, "   inner_box: ", pfa_knots.inner_box)
    println(io, "   quadratic: ", pfa_knots.quadratic)
    println(io, "   distortion: ", pfa_knots.distortion)
    println(io, "   image_plane: ", pfa_knots.image_plane)
    println(io, "   enu_align: ", pfa_knots.enu_align)
    return nothing
end
