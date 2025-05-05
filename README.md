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

After installing, verify that OceanLight works as intended by:

```Julia
Pkg.test("OceanLight")
```

## Result

<img  src="https://raw.githubusercontent.com/haoboatlab/OceanLight.jl/main/docs/src/assets/Center1e8.png" width="600" align="center">
*Simulation of 1e8 Photons at the center of the flat surface. (a) irradiance field on the horizontal plane at 30 m depth. (b) irradiance field on the horizontal plane at 150 m depth. (c) irradiance field on the vertical plane at the center.*


<img  src="https://raw.githubusercontent.com/haoboatlab/OceanLight.jl/main/docs/src/assets/Wholegrid1000.png" width="600" align="center">
*Simulation of 1000 Photons at the every grid point with observed surface elevation. (a) irradiance field on the horizontal plane at 30 m depth. (b) irradiance field on the horizontal plane at 150 m depth. (c) irradiance field on the vertical plane at the center.*


## Reference 

Hao, X., & Shen, L. (2022). A novel machine learning method for accelerated modeling of the downwelling irradiance field in the upper ocean. Geophysical Research Letters, 49, e2022GL097769. https://doi.org/10.1029/2022GL097769

Kirk, J. T. O. (1981). Monte Carlo procedure for simulating the penetration of light into natural waters. In Technical paper - Commonwealth Scientific and Industrial Research Organization (Vol. 36).
