# OceanLight.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://haoboatlab.github.io/OceanLight.jl/dev/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://mit-license.org/)

## Overview

**OceanLight.jl** calculates the downwelling irradiance field in the upper ocean. By implementing the Monte Carlo method, this simulation achieves the path of Photons: starting from the refraction in the air-water interface, to some specific depth underneath the water body. 

## Installation 

OceanLight requires Julia software. To do so, 

1. [Install Julia](https://julialang.org/downloads/) 

2. Lauch Julia and type 

```Julia
using Pkg
Pkg.add("OceanLight")
```

## Result

<img  src="https://raw.githubusercontent.com/haoboatlab/OceanLight.jl/main/docs/src/assets/center1e7.png" width="600" align="center">
*Simulation of 1e7 Photons at the center with no surface elevation, in 80 m by 80 m by 190 m physical size with 512 by 512 by 190 resolution, absorbtance coefficent a = 0.0196 and scattering coefficient b = 0.0031. Sub-Figure 1: X-Y cross section at 30 m depth. Sub-Figure 2: X-Y cross section at 150 m depth. Sub-Figure 3: X-Z cross section at the center.*


<img  src="https://raw.githubusercontent.com/haoboatlab/OceanLight.jl/main/docs/src/assets/WholeGrid100.png" width="600" align="center">
*Simulation of 100 Photons at the every grid point with observed surface elevation, in 80 m by 80 m by 190 m physical size with 512 by 512 by 190 resolution, absorbtance coefficent a = 0.0196 and scattering coefficient b = 0.0031. Sub-Figure 1: X-Y cross section at 30 m depth. Sub-Figure 2: X-Y cross section at 150 m depth. Sub-Figure 3: X-Z cross section at the center.*


## Reference 

Kirk, J. T. O. (1981). Monte Carlo procedure for simulating the penetration of light into natural waters. In Technical paper - Commonwealth Scientific and Industrial Research Organization (Vol. 36).
