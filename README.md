# RadarData

[![Documentation](https://img.shields.io/badge/docs-public-blue.svg)](https://matrix-research-inc.github.io/RadarData.jl/dev/)
[![Documentation](https://img.shields.io/badge/docs-internal-orange.svg)](https://pages-git.matrixresearch.com/programs/alg/julia/radardata)
[![Latest Release](https://img.shields.io/github/v/release/matrix-research-inc/RadarData.jl)](https://github.com/matrix-research-inc/RadarData.jl/releases/tag/public-release-v0.7.3)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A [Julia](https://julialang.org/) package for handling synthetic aperture radar (SAR) data types. 

## Features

- **Video Phase History (VPH)**: Core data structure for range-compressed radar pulse data
- **SAR Image Format**: Georeferenced SAR image container with coordinate transformations
- **Position, Navigation, Timing (PNT)**: Platform state management and interpolation
- **Coordinate Transformations**: Support for ECEF, ENU, LLA, and horizon coordinate systems
- **Range Compression**: Chirp synthesis and pulse compression utilities
- **Polar Format Algorithm**: PFA metadata structures for SAR image formation
- **Package Extensions**: Optional conversion support for other formats

## License

This project is licensed under the **GNU Affero General Public License v3.0 (AGPLv3)**.

You are free to use, modify, and distribute this software under the terms of the AGPLv3.
If you run a modified version of this software as a network service, you must make the source code of your modified version available to users.

**Commercial License:** Matrix Research, Inc. offers this software under a separate commercial license for organizations that wish to use it in proprietary applications without the AGPLv3 obligations. 
Contact us for commercial licensing options.

See the [LICENSE](LICENSE) file for full license text.

## Contributing

We welcome contributions! 
Before submitting a pull request, please note:

- **Contributor License Agreement (CLA)**: External contributors must sign a CLA granting Matrix Research, Inc. the right to use and sublicense your contributions. 
This allows us to offer both open-source and commercial licenses.
- **Code Review**: All contributions undergo review for compliance with our intellectual property and security policies.

To get started or for questions about contributing, contact Sam Pine <sam.pine@matrixresearch.com>.

## Support

Please contact Sam Pine <sam.pine@matrixresearch.com> for support.

## Acknowledgements
The authors thank AFRL/Sensors Directorate for their support of this work via SBIR Contract FA2377-24-C-B023.
This software is approved for public release. 
See the [public release history](docs/public_release_history.md) for details.