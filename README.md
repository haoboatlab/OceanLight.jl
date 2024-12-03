# OceanLight.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://haoboatlab.github.io/OceanLight.jl/dev/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://mit-license.org/)

## Overview

**OceanLight.jl** calculates the downwelling irradiance field in the upper ocean. By implementing the Monte Carlo method, this simulation achieves the path of Photons: starting from the refraction in the air-water interface, to some specific depth underneath the water body. 

## Installation 

OceanLight requires Julia software. To do so, 

1. [Download Julia](https://julialang.org/downloads/) 

2. Lauch Julia and type 

```@example
using Pkg
Pkg.add("OceanLight")
```

## Reference 

Kirk, J. T. O. (1981). Monte Carlo procedure for simulating the penetration of light into natural waters. In Technical paper - Commonwealth Scientific and Industrial Research Organization (Vol. 36).
