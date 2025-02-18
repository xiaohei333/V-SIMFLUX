# V-SIMFLUX

<!-- This should be updated after the paper is published. -->
<!-- [![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.optlaseng.2025.108851-blue.svg)](https://doi.org/10.1016/j.optlaseng.2025.108851) -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

This repository contains the official implementation of the research described in:

> Enhancing precision for simultaneous 3D localization and 3D orientation with structured illumination
>
> [![Zhengyi Zhan](https://img.shields.io/badge/Zhengyi%20Zhan-181717?logo=github&logoColor=white)](https://github.com/ZhengyiZ) [![Xiaofan Sun](https://img.shields.io/badge/Xiaofan%20Sun-181717?logo=github&logoColor=white)](https://github.com/xiaohei333) [![Chenying He](https://img.shields.io/badge/Chenying%20He-181717?logo=github&logoColor=white)](https://github.com/Haibara647) et al.
<!-- This should be updated after username changed. -->

*Any reuse of this code should cite the original associated publication.*

## Introduction

V-SIMFLUX is an advanced single-molecule orientation localization microscopy (SMOLM) method that uniquely combines sequential structured illumination (SSI) with the Vortex point spread function (PSF). This integration typically enables:

- 2.9× improvement in lateral localization precision
- 1.7× enhancement in azimuthal orientation determination

The method provides a comprehensive framework for simultaneous 3D position and orientation determination of single molecules, particularly effective in limited-photon scenarios.

## Quick Start

### Dataset Visualization

Navigate to the `Dataset` folder and run the corresponding figure script (e.g. `Fig_1_Demo/DispA.m`) to visualize specific results:

- **Fig_1_Demo**: Demonstration results
- **Fig_2_S4_S5_S6_CRLB**: CRLB analysis
- **Fig_3_SBR**: Signal-to-background ratio analysis
- **Fig_4_Ring**: Ring structure analysis
- **Fig_S2_PSF, Fig_S3_NoiseImage**: Supplementary data

### Running the Code

The code is organized in the following structure:

- **funextra**: Utility functions for image loading and peak finding
- **funfit**: Core fitting and analysis functions
- **funnat**: Nodal Aberration Theory related functions
- **V-SIMFLUX**: V-SIMFLUX specific implementations

<!-- Simulation workflow:

1. Configure parameters:
   - `set_parameters_vortex_sim.m` - For simulated data
   - `set_parameters_vortex_lambdaDNA.m` - For DNA samples
   - `set_parameters_zstack_bead.m` - For bead calibration

2. Run analysis:
   - `vecfitcpu_vortex_simfits.m` - Simulation analysis
   - `vecfitcpu_vortex_lambdaDNA.m` - DNA sample analysis
   - `vecfitcpu_zstack_bead.m` - Bead calibration -->

## Citation & Reference

<!-- This should be updated after the paper is published. -->
<!-- If you use this work in your research, please cite:

``` bibtex
xxx
``` -->

This project builds upon *Vortex PSF* [![DOI](https://img.shields.io/badge/DOI-10.1038/s41467--021--26228--5-blue)](https://doi.org/10.1038/s41467-021-26228-5)

``` bibtex
@article{hullemanSimultaneousOrientation3D2021,
    author = {Hulleman, Christiaan N. and Thorsen, Rasmus {\O} and Kim, Eugene and Dekker, Cees and Stallinga, Sjoerd and Rieger, Bernd},
    title = {Simultaneous Orientation and 3D Localization Microscopy with a Vortex Point Spread Function},
    year = {2021},
    doi = {10.1038/s41467-021-26228-5},
    journal = {Nature Communications}
}
```

## License

This project incorporates code from [![Vortex PSF](https://img.shields.io/badge/vecfitcpu__vortex-0076A8?logo=github)](https://github.com/imphys/vecfitcpu_vortex), which is licensed under the [Apache-2.0 License](./Code/LICENSE).

V-SIMFLUX modifications are provided under the [MIT License](./Code/V-SIMFLUX/LICENSE).
