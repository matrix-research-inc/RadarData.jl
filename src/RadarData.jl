#=
    License ID: SEAL_B
=#
module RadarData

using StaticArrays: SVector, SMatrix, MMatrix
using StatsBase: mean
using FFTW: fftshift, fft, fft!, ifftshift, ifft, ifft!
using DSP: unwrap
using Geodesy: ECEF, ENU, LLA, wgs84, wgs84_ellipsoid
using Interpolations: interpolate, extrapolate, Gridded, Linear, NoInterp
using Rotations: RotXYZ, RotZX, QuatRotation
using LinearAlgebra
using Printf
using DocStringExtensions
using KernelAbstractions
using KernelHelper

"""
    c0 = 299792458

Speed of light in a vaccuum, expressed in meters per second
"""
const c0 = 299792458
export c0   # Speed of light (m/s)

include("fft_shortcuts.jl")
include("coordinates.jl")
include("pnt.jl")
include("vph.jl")
include("test_data.jl")
include("range_compress.jl")
include("sar_image.jl")
include("image_metadata.jl")

end # module RadarData
