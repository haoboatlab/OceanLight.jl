# OceanLight.jl Documentation

## Overview 

**OceanLight.jl** calculates the downwelling irradiance field in the upper ocean. It uses the Monte Carlo method to simulate the trajectories of photons, from their refraction at the air–water interface to specific depths beneath the surface [^1] [^2]. 

Optical oceanography concerns all aspects of light and its interaction with seawater, which are crucial for addressing problems related to physical, biological, and chemical oceanographic processes, such as phytoplankton photosynthesis, biogeochemical cycles, and climate change [^3] [^4]. However, due to the complex interaction between light and free-surface wave geometry, the irradiance distribution can be highly variable [^5] [^6], making it difficult to obtain analytical solutions. **OceanLight.jl** enables physics-based, reproducible simulations of light distribution beneath the ocean surface.

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

![Center1e8](https://raw.githubusercontent.com/haoboatlab/OceanLight.jl/main/docs/src/assets/Center1e8.png)

*Simulation of $10^{8}$ Photons at the center of the flat surface. (a) irradiance field on the horizontal plane at $30\ \mathrm{m}$ depth. (b) irradiance field on the horizontal plane at $150\ \mathrm{m}$ depth. (c) irradiance field on the vertical plane at the center*

![Wholegrid1000](https://raw.githubusercontent.com/haoboatlab/OceanLight.jl/main/docs/src/assets/Wholegrid1000.png)

*Simulation of 1000 Photons at the every grid point with observed surface elevation. (a) irradiance field on the horizontal plane at $30\ \mathrm{m}$ depth. (b) irradiance field on the horizontal plane at $150\ \mathrm{m}$ depth. (c) irradiance field on the vertical plane at the center.*

## Reference 

[^1]: Hao, X., & Shen, L. (2022). A novel machine learning method for accelerated modeling of the downwelling irradiance field in the upper ocean. Geophysical Research Letters, 49, e2022GL097769. https://doi.org/10.1029/2022GL097769

[^2]: Kirk, J. T. O. (1981). Monte Carlo procedure for simulating the penetration of light into natural waters. In Technical paper - Commonwealth Scientific and Industrial Research Organization (Vol. 36).

[^3]: Dickey, T. D., Kattawar, G. W., & Voss, K. J. (2011). Shedding new light on light in the ocean. Physics Today, 64(4), 44-49.

[^4]: Dickey, T., Lewis, M., & Chang, G. (2006). Optical oceanography: recent advances and future directions using global remote sensing and in situ observations. Reviews of geophysics, 44(1).

[^5]: Darecki, M., Stramski, D., & Sokólski, M. (2011). Measurements of high‐frequency light fluctuations induced by sea surface waves with an Underwater Porcupine Radiometer System. Journal of Geophysical Research: Oceans, 116(C7).

[^6]: Gernez, P., Stramski, D., & Darecki, M. (2011). Vertical changes in the probability distribution of downward irradiance within the near‐surface ocean under sunny conditions. Journal of Geophysical Research: Oceans, 116(C7).