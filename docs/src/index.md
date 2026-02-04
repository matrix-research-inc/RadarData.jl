```@meta
CurrentModule = RadarData
```

# RadarData

[RadarData](https://git.matrixresearch.com/Programs/alg/julia/RadarData.jl) provides data structures and utilities for working with synthetic aperture radar (SAR) data. 
The package handles radar phase history data, position/navigation/timing (PNT) information, coordinate transformations, and motion compensation.

## Overview

The main components of RadarData.jl are:

- **VPH (Video Phase History)**: Core data structure for radar returns
- **PNT (Position, Navigation, Timing)**: Platform state information
- **Coordinate Systems**: Transformations between ECEF, ENU, LLA, and horizon coordinates
- **Motion Compensation**: Phase correction for platform motion
- **Range Compression**: Pulse compression utilities
- **Test Data Generation**: Synthetic data for algorithm development

## Constants

```@docs
c0
```

## Video Phase History (VPH)

The core data type maintained by RadarData.jl is the video phase history, or VPH. 
This structure contains range-compressed radar returns along with associated metadata including platform positions, timing, and frequency information.

### VPH Type

```@docs
VPH
VPH(::AbstractMatrix{Complex{T}} where T<:AbstractFloat,
    ::AbstractVector,
    ::AbstractVector,
    ::AbstractVector,
    ::AbstractMatrix)
VPH(::Int, ::Int)
```

### VPH Metadata Functions

These functions extract or compute metadata from VPH objects:

```@docs
pri(::VPH)
get_freq_list
spatial_freqs(::VPH)
bandwidth(::VPH)
center_freq(::VPH)
range_res(::VPH)
velocity_res(::VPH)
range_axis(::VPH)
fast_time_axis(::VPH)
velocity_axis(::VPH)
time_axis(::VPH)
doppler_axis(::VPH)
cross_range_extent(::VPH)
bistatic_angles(::VPH)
bisectors(::VPH)
range_dir_enu(::VPH)
cross_range_dir_enu(::VPH)
bp_layover_projection(::VPH)
azimuth_rad(::VPH)
elevation_rad(::VPH)
horizon_coords(::VPH)
show
```

### VPH Domain Conversion

Convert between frequency domain and range (fast-time) domain:

```@docs
vph_to_rs!(::VPH)
vph_to_rs(::VPH)
rs_to_vph!(::VPH)
rs_to_vph(::VPH)
```

### VPH Data Selection

Extract subsets of VPH data:

```@docs
freq_truncate!(::VPH,::UnitRange{Int})
freq_truncate!(::VPH,::Real,::Real)
time_truncate!(::VPH,::AbstractVector{Int})
time_truncate!(::VPH,::Real,::Real)
extract_cpi(::VPH, ::AbstractVector{Int})
```

## Position, Navigation, and Timing (PNT)

PNT structures store platform state information including position, velocity, and attitude.

### PNT Types

```@docs
PVAHeader
PVA
PNT
PNTHistory
```

### PNT Functions

```@docs
receive_pva!
interp_states
interp_ecef
range_to_point
range_rate_to_point
```

## Coordinate Transformations

Functions for working with coordinate systems. 
VPH objects are always expressed in a local coordinate system, assumed to 
be the East-North-Up (ENU) tangent plane to the WGS84 ellipsoid at `vph.lla`.

```@docs
cart_to_horizon
cart_to_horizon_unwrap
horizon_to_cart
azimuth_rad
elevation_rad
slant_plane_rot_2D
slant_plane_rot_3D
ned_to_enu
body_to_ecef
enu_rotation
intersect_los_plane
```

## Motion Compensation

Correct for platform motion effects:

```@docs
phase_to_range
range_to_phase
motion_compensate!
motion_compensate
```

## Range Compression

Pulse compression utilities:

```@docs
synth_chirp
synth_chirp!
range_compress
range_compress!
```

## Test Data Generation

Generate synthetic radar data for testing:

```@docs
ScatterCenters
generate_vph
```

## FFT Utilities

Convenience functions for FFT operations:

```@docs
fft_spacing
sfft
sifft
vph_to_rd
rd_to_vph
```

## Support for Conventions

Common conventions are encapsulated in function handles for easy recall:

```@docs
center_freq
spatial_freqs
bandwidth
range_res
velocity_res
range_axis
velocity_axis 
```

## Index
```@index
```
