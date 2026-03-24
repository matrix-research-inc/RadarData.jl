
export SARImage,
       resolution,
       srp_ecef,
       image_corners_ecef,
       image_corners_lla,
       enu_frame,
       image_axes

"""
SARImage

Structure for containing SAR image data with axes and labels.

$(TYPEDFIELDS)
"""
@kwdef mutable struct SARImage{T <: AbstractFloat}
    "Image pixel data. Expected to be either ComplexF32 or ComplexF64"
    data::Matrix{Complex{T}}
    "Location of pixel grid centers, increasing down first dimension (row)"
    row_axis::Vector{Float64}
    "Location of pixel grid centers, increasing across second dimension (column)"
    col_axis::Vector{Float64}
    "Label for row axis"
    row_label::String
    "Label for column axis"
    col_label::String
    "Scene reference point in LLA coordinates"
    srp_lla::SVector{3, Float64}
    "Pixel corresponding to scene reference point"
    srp_pixel::SVector{2, Int} = collect(fld.(size(data), 2) .+ 1)
    "Unit vector pointing in direction of increasing first direction in ECEF coordinate frame"
    row_unit_vector::SVector{3, Float64}
    "Unit vector pointing in direction of increasing second direction in ECEF coordinate frame"
    col_unit_vector::SVector{3, Float64}
    "Identifier for algorithm used to form the image"
    algorithm::String
    "Full image size, if this image is a sub-image"
    full_size::SVector{2, Int64} = SVector{2}(size(data, 1), size(data, 2))
    "Row and column index of the first pixel in the sub-image, relative to the full image"
    start_idxs::SVector{2, Int64} = SVector{2}(1, 1)
end

Base.size(image::SARImage) = size(image.data)
Base.size(image::SARImage, dims::Int) = size(image.data, dims)
Base.eltype(image::SARImage) = eltype(image.data)
Base.length(image::SARImage) = length(image.data)

"""
    resolution(image)

Return static vector with pixel spacing in meters of SARImage in row and column axes.
"""
resolution(image::SARImage) = SVector{2, Float64}(
    image.row_axis[2] - image.row_axis[1],
    image.col_axis[2] - image.col_axis[1]
)

function Base.getindex(image::SARImage, rows, cols)
    return SARImage(;
        data = image.data[rows, cols],
        row_axis = image.row_axis[rows],
        col_axis = image.col_axis[cols],
        row_label = image.row_label,
        col_label = image.col_label,
        srp_lla = image.srp_lla,
        srp_pixel = image.srp_pixel,
        row_unit_vector = image.row_unit_vector,
        col_unit_vector = image.col_unit_vector,
        algorithm = image.algorithm,
        full_size = image.full_size,
        start_idxs = image.start_idxs .+ SVector{2, Int64}(rows[1] - 1, cols[1] - 1)
    )
end

function srp_ecef(image::SARImage)
    ECEF(LLA(image.srp_lla[1], image.srp_lla[2], image.srp_lla[3]), wgs84)
end

"""
    image_axes(image)

Return image axes, expressed in relative ECEF coordinates.
Column order is down range, increasing cross range, and up out of imaging plane.
"""
function image_axes(image::SARImage)
    ecef_axes = zeros(MMatrix{3, 3, Float64})
    ecef_axes[:, 1] = image.row_unit_vector
    ecef_axes[:, 2] = image.col_unit_vector
    ecef_axes[:, 3] = cross(image.row_unit_vector, image.col_unit_vector)
    return ecef_axes
end

"""
    enu_frame(image)

Return East-North-Up (ENU) coordinate frame centered on `image` SRP, expressed in relative ECEF coordinates.
"""
function enu_frame(image::SARImage)
    srp = srp_ecef(image)
    enu_axes = zeros(MMatrix{3, 3, Float64})
    enu_axes[:, 1] = ECEF(ENU(1, 0, 0), srp, wgs84) - srp
    enu_axes[:, 2] = ECEF(ENU(0, 1, 0), srp, wgs84) - srp
    enu_axes[:, 3] = ECEF(ENU(0, 0, 1), srp, wgs84) - srp
    return enu_axes
end

"""
    image_corners_ecef(image)

Return corners of image in ECEF coordinates.
Order of columns is first row, first column; first row, last column; last row, last column; last row, first column.
Note, this does not account for Earth curvature.
"""
function image_corners_ecef(image::SARImage)
    center = srp_ecef(image)
    corners = zeros(3, 4)
    corners[:, 1] = center .+ image.row_axis[1] .* image.row_unit_vector .+
                    image.col_axis[1] .* image.col_unit_vector
    corners[:, 2] = center .+ image.row_axis[1] .* image.row_unit_vector .+
                    image.col_axis[end] .* image.col_unit_vector
    corners[:, 3] = center .+ image.row_axis[end] .* image.row_unit_vector .+
                    image.col_axis[end] .* image.col_unit_vector
    corners[:, 4] = center .+ image.row_axis[end] .* image.row_unit_vector .+
                    image.col_axis[1] .* image.col_unit_vector
    return corners
end

"""
    image_corners_lla(image)

Return corners of image in LLA coordinates.
Order of columns is first row, first column; first row, last column; last row, last column; last row, first column.
Note, this does not account for Earth curvature.
"""
function image_corners_lla(image::SARImage)
    ecef_corners = image_corners_ecef(image)
    lla_corners = zeros(3, 4)
    for idx in axes(ecef_corners, 2)
        lla_coord = LLA(ECEF(ecef_corners[:, idx]), wgs84)
        lla_corners[1, idx] = lla_coord.lat
        lla_corners[2, idx] = lla_coord.lon
        lla_corners[3, idx] = lla_coord.alt
    end
    return lla_corners
end

"""
	show(io, image)

Custom display method for SARImage.
"""
function Base.show(io::IO, image::SARImage{T}) where {T <: AbstractFloat}
    println(io, size(image, 1), "x", size(image, 2), " SARImage{", T, "}:")
    println(
        io,
        "    data: ",
        size(image.data),
        @sprintf(" [%.2f%+.2fim %.2f%+.2fim ... %.2f%+.2fim %.2f%+.2fim]",
            real(image.data[1]),
            imag(image.data[1]),
            real(image.data[2]),
            imag(image.data[2]),
            real(image.data[end - 1]),
            imag(image.data[end - 1]),
            real(image.data[end]),
            imag(image.data[end]))
    )
    println(
        io,
        "    row_axis: ",
        size(image.row_axis),
        @sprintf(" [%.2f ... %.2f]", image.row_axis[1], image.row_axis[end])
    )
    println(
        io,
        "    col_axis: ",
        size(image.col_axis),
        @sprintf(" [%.2f ... %.2f]", image.col_axis[1], image.col_axis[end])
    )
    println(io, "    row_label: ", image.row_label)
    println(io, "    col_label: ", image.col_label)
    println(
        io,
        "    srp_lla: ",
        @sprintf("[%.2f %.2f %.2f]", image.srp_lla[1], image.srp_lla[2], image.srp_lla[3])
    )
    println(
        io,
        "    srp_pixel: ",
        @sprintf("[%.0f %.0f]", image.srp_pixel[1], image.srp_pixel[2])
    )
    println(
        io,
        "    row_unit_vector: ",
        @sprintf("[%.2f %.2f %.2f]",
            image.row_unit_vector[1],
            image.row_unit_vector[2],
            image.row_unit_vector[3])
    )
    println(
        io,
        "    col_unit_vector: ",
        @sprintf("[%.2f %.2f %.2f]",
            image.col_unit_vector[1],
            image.col_unit_vector[2],
            image.col_unit_vector[3])
    )
    println(io, "    algorithm: ", image.algorithm)
    println(
        io,
        "    full_size: ",
        @sprintf("[%.0f %.0f]", image.full_size[1], image.full_size[2])
    )
    println(
        io,
        "    start_idxs: ",
        @sprintf("[%.0f %.0f]", image.start_idxs[1], image.start_idxs[2])
    )
end
