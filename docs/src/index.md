# OceanLight.jl Documentation

## Overview 

**OceanLight.jl** calculates the downwelling irradiance field in the upper ocean. By implementing the Monte Carlo method, this simulation achieves the path of Photons: starting from the refraction in the air-water interface, to some specific depth underneath the water body. 

## Quick Install 

1. [Install Julia](https://julialang.org/downloads/)

2. Launch Julia and type

```julia
using Pkg
Pkg.add("OceanLight")
```

After installing, verify that OceanLight works as intended by:

```Julia
Pkg.test("OceanLight")
```

## Result

<img  src="assets/Center1e8.png" width="600" align="center">

*Simulation of $10^{8}$ Photons at the center of the flat surface. (a) irradiance field on the horizontal plane at $30\ \mathrm{m}$ depth. (b) irradiance field on the horizontal plane at $150\ \mathrm{m}$ depth. (c) irradiance field on the vertical plane at the center.*

<img  src="assets/Wholegrid1000.png" width="600" align="center">

*Simulation of 1000 Photons at the every grid point with observed surface elevation. (a) irradiance field on the horizontal plane at $30\ \mathrm{m}$ depth. (b) irradiance field on the horizontal plane at $150\ \mathrm{m}$ depth. (c) irradiance field on the vertical plane at the center.*