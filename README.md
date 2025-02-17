# V-SIMFLUX

<!-- This should be updated after the paper is published. -->
<!-- This repository contains the implementation of V-SIMFLUX, accompanying the paper [*Enhancing precision for simultaneous 3D localization and 3D orientation with structured illumination*]() by Xiaofan Sun, Zhengyi Zhan, Chenying He, et al. -->

<!-- This should be updated after the paper is published. -->
<!-- [![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.optlaseng.2025.108851-blue.svg)](https://doi.org/10.1016/j.optlaseng.2025.108851) -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

Any reuse of this code should cite the original associated publication.

## Introduction

V-SIMFLUX is a novel method that combines sequential structured illumination (SSI) with the Vortex point spread function (PSF) for enhanced single-molecule localization microscopy. By integrating these approaches, V-SIMFLUX achieves:

- 2.9× improvement in lateral localization precision
- 1.7× enhancement in azimuthal orientation determination

This implementation provides a complete framework for simultaneous 3D position and orientation determination of single molecules, addressing the key challenges in limited-photon scenarios.

## Quick Start

### Visualizing Data

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

This project builds upon two published methods:

- **Vortex method** [![DOI](https://img.shields.io/badge/DOI-10.1038/s41467--021--26228--5-blue)](https://doi.org/10.1038/s41467-021-26228-5)

    ``` bibtex
    Hulleman, Christiaan N., et al. "Simultaneous orientation and 3D localization microscopy with a Vortex point spread function." Nature Communications 12.1 (2021): 5934.
    ```

- **SIMFLUX method** [![DOI](https://img.shields.io/badge/DOI-10.1038/s41592--019--0657--7-blue)](https://doi.org/10.1038/s41592-019-0657-7)

    ``` bibtex
    Cnossen, Jelmer, et al. "Localization microscopy at doubled precision with patterned illumination." Nature methods 17.1 (2020): 59-63.
    ```

## License

This project incorporates code from the following sources:

- [vecfitcpu_vortex](https://github.com/imphys/vecfitcpu_vortex)
  
  Original Vortex implementation - Licensed under the [Apache-2.0 License](./Code/LICENSE).
  
- [simflux](https://github.com/qnano/simflux)
  
  Original SIMFLUX implementation - Licensed under the [MIT License](https://github.com/qnano/simflux?tab=MIT-1-ov-file#readme).

Our V-SIMFLUX modifications are licensed under the [MIT License](./LICENSE).
