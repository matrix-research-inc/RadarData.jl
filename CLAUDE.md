# CLAUDE.md
This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview
RadarData is a Julia package for handling synthetic aperture radar (SAR) data types. 
The core data structure is Video Phase History (VPH), which represents range-compressed radar pulse data with associated position, navigation, and timing (PNT) metadata.

## Architecture

### Core Data Structures

**VPH (Video Phase History)**: The central data structure (`src/vph.jl`)
- `data`: Complex-valued matrix (frequency or fast-time × slow-time)
- `freq_list_hz`: Frequency list in Hz
- `ref_range_m`: Two-way range to center range bin per pulse
- `slow_time_s`: Per-pulse collection time
- `tx_pos_enu`, `rx_pos_enu`: Transmitter/receiver positions (East-North-Up coordinates)
- `srp_enu`, `srp_lla`: Scene reference point in ENU and LLA coordinates
- `range_domain`: Boolean flag indicating if data is in frequency or range (fast-time) domain
- Supports both ComplexF32 and ComplexF64 data, but position/timing metadata is always Float64

**PNT (Position, Navigation, Timing)**: Platform state information (`src/pnt.jl`)
- `PVA`: Position, Velocity, Attitude message following ASPN convention
- `PNTHistory`: Collection of PVA messages over time
- Functions for interpolating states, computing ranges, and coordinate transformations

### Module Organization

The package is organized into focused modules included from `src/RadarData.jl`:
- `coordinates.jl`: Coordinate system conversions (Cartesian ↔ horizon, azimuth/elevation, slant plane geometry)
- `pnt.jl`: Position, navigation, timing structures and interpolation
- `vph.jl`: Video Phase History structure and domain transformations
- `motion_compensate.jl`: Phase/range conversions and motion compensation algorithms
- `range_compress.jl`: Chirp synthesis and range compression operations
- `fft_shortcuts.jl`: FFT utilities for radar processing
- `test_data.jl`: Synthetic data generation for testing

### Package Extensions

RadarData uses Julia's package extension system for optional functionality (requires Julia 1.9+):

1. **MTR** (`ext/MTR.jl`): DARPA Moving Target Recognition legacy data format support
   - Activated when `MAT` package is loaded
   - Provides `read_vph_legacy()` for MATLAB-compatible VPH files

2. **TSC** (`ext/TSC.jl`): Time-space-coherent data format
   - Activated when `MAT` package is loaded
   - Functions for reading TSC-format radar data

3. **SoulExt** (`ext/SoulExt/`): RFSoC interface integration
   - Activated when `SoulHardware`, `SoulINS`, and `BEVE` packages are loaded
   - Converts RxData types to VPH format
   - Handles BEVE file deserialization from RFSoC hardware

### Coordinate Systems

The package uses multiple coordinate systems:
- **ECEF**: Earth-Centered, Earth-Fixed (global Cartesian)
- **LLA**: Latitude, Longitude, Altitude (geodetic)
- **ENU**: East-North-Up (local Cartesian, default for VPH positions)
- **NED**: North-East-Down (used for body frame attitude)
- **Horizon**: Range, azimuth, elevation (spherical)

Conversion functions are in `src/coordinates.jl` and `src/pnt.jl`.

### Key Conventions

- **In-place operations**: Functions ending with `!` modify their first argument
- **Dual domain support**: VPH data can be in frequency or range (fast-time) domain
  - Use `vph_to_rs!()` / `rs_to_vph!()` to transform between domains
  - Check `vph.range_domain` flag to determine current domain
- **Speed of light**: Effective speed `c_eff_ms` accounts for atmospheric propagation
  - Default vacuum speed available as constant `c0 = 299792458` m/s
- **FFT convention**: Controlled by `kernel_sign` field (default -1)
- **License headers**: All source files include `#= License ID: SEAL_B =#` header

## Development Notes

- The package targets Julia 1.12.4+ (see `Project.toml`)

