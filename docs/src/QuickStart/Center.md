# Quick Start 

In this example, we introduce the `OceanLight` calculation of downwelling irradiance field underneath the flat surface. The code example here can be directly pasted onto Julia terminal, run through `.jl` script file, or IJulia notebook. 

First, we import `OceanLight` packages and another `random` dependent package, that will be used as a seed to generate random number in the Monte Carlo simulation.

```@setup  Center
cd(mktempdir()) 
using Pkg 
Pkg.add("OceanLight") 
Pkg.add("Plots") 
```
```@example Center 
using OceanLight 
using Random
```

## Problem

In this example, the problem is to calculate the downwelling irradiance field, when the surface is completely flat, and a total of 10,000,000 photons is focused at a single point in the center. Our domain of interest is defined as $x,y \in \left[\mathrm{-10m},\mathrm{10m}\right]$, and $z \in \left[\mathrm{-190m},\mathrm{10m}\right]$ in depth, corresponding to a grid resolution of $512 \times 512 \times 200$  points. Periodic boundary conditions are applied at the domain boundaries. The attenuation properties of water are characterized by an absorption coefficient of $a = 0.0196$ and scattering coefficient $b = 0.0031$, which corresponding to the optical properties of sea water at wavelength $490 \mathrm{nm}$ [^1]. 

##  Initial Condition 

OceanLight accesses all input variables in `.yml` format and stores their values in `Param` structure.

All input variables required by OceanLight can be separated into 3 categories:
1. **Irradiance:** Setting up the dimension of downwelling irradiance field calculation grid. 
2. **Photon:** Number of Photons and optical properties of the water.
3. **Wave:** Structure of surface wave field. 

```@example Center
# irradiance
nz = 200                    # Number of total grid point in z direction
dz = 1                      # Physical length of grid spacing in z direction
nxe = 512                   # Number of grid point of the calculation grid in x direction  
nye = 512                   # Number of grid point of the calculation grid in y direction  
num = 31                    # Constant value (number of angle measurement in Kirk,1981)
ztop = 10                   # Number of grid point in air phase in z direction
# photon
nphoton = 10000000          # Number of photons being generated at each grid point 
kr = 10                     # Dummy variable (in developement, not being used)                    
nxp = 512                   # Number of grid points in x direction where photon can be emitted
kbc = 0                     # Binary value 0 and 1 depending on Boundary condition being implemented 
b = 0.0031                  # Scattering coefficient 
nyp = 512                   # Number of grid points in x direction where photon can be emitted
a = 0.0196                  # Absorbtance coefficient
# wave
pey = 2*pi/20.0             # Lowest wavenumber that can be captured during the derivative of surface elevation in x direction
nxeta = 512                 # Number of surface elevation grid point in x direction
nyeta = 512                 # Number of surface elevation grid point in y direction
pex = 2*pi/20.0             # Lowest wavenumber that can be captured during the derivative of surface elevation in y direction
```
**NOTE:** `num` must be set to a constant value of 31 (number of angle measurement in Kirk, 1981 [^2]). `kbc` accepts only binary values (0 or 1), representing periodic boundary conditions. In contrast to grid spacing `dz` in z-direction, the grid spacings in the x- and y-directions (`dx` and `dy`) are calculated automatically by `OceanLight.readparams`,  using the formulas $dx \times nxe  = \frac{2\pi}{pex}$, and $dy \times nye = \frac{2\pi}{pey}$. For more detail on all parameters  used, see [`Simulation parameters`](@ref simulation_parameters). 

To create a input variable file suitable for this package, user can either create a new file in `.yml` format, copy, and paste the code block above, or using a build-in function `OceanLight.writeparams` to automate the process. 

The function `OceanLight.writeparams` converts the dictionary of input variables into the `light.yml` file. 

```@example Center
data=Dict("irradiance"=>Dict("nxe"=>nxe,"nye"=>nye,"nz"=>nz,"dz"=>dz,"ztop"=>ztop,"num"=>num),
            "wave"=>Dict("pex"=>pex,"pey"=>pey,"nxeta"=>nxeta,"nyeta"=>nyeta),
            "photon"=>Dict("nxp"=>nxp,"nyp"=>nyp,"nphoton"=>nphoton,"a"=>a,"b"=>b,"kr"=>kr,"kbc"=>kbc))

OceanLight.writeparams(data)
```

Once the input variable file is generated, this package calculates other related variables and stores all values in the `Param` structure, through `OceanLight.readparams`. The default input variable file name that will be read is `light.yml`, but can be specified by users. In this case, we store `Param` structure in variable `p`. 

```@example Center
p = OceanLight.readparams()
```

## Initialize Parameters

Before the simulation, user needs to declare and initialize all parameters and their dimensions. 

During the light refraction between two mediums (air and water) calculation, `OceanLight` requires the surface elevation attribution: $\eta$; surface elevation, $\eta_{x}$; partial derivative of $\eta$ in x direction, and $\eta_{y}$; partial derivative of $\eta$ in y direction. All surface elevation distribution, as specified above, need to have the same dimension with the incoming photons' grid. Hence, 

```@example Center
η = zeros(p.nxs,p.nys)
ηx = zeros(p.nxs,p.nys)
ηy = zeros(p.nxs,p.nys)
```
After photons' interaction with the surface, OceanLight requires the information of specific coordinate of photon in cartesian grid ${xpb,ypb,zpb}$, the direction in which photon will travel in polar coordinate ${θ,ϕ}$, and the fraction of light that transmit through the water ${fres}$: all in the dimension of incoming photon grid size. 

```@example Center
xpb = zeros(p.nxp,p.nyp);
ypb = zeros(p.nxp,p.nyp);
zpb = zeros(p.nxp,p.nyp);
θ = zeros(p.nxp,p.nyp);
ϕ = zeros(p.nxp,p.nyp);
fres = zeros(p.nxp,p.nyp);
```

During the air-water interaction process, OceanLight simulates the photons transfer directly downward from the air side, interacts with the water surface, and transfer down into water medium. 

User can generate random surface elevation attribution $\{\eta,\eta_{x},\eta_{y}\}$ with `OceanLight.setwave!`, or provided specific data $\{\eta_{0},\eta_{x0},\eta_{y0}\}$  . OceanLight can map the user's provided data of $\{\eta_{0},\eta_{x0},\eta_{y0}\}$, which might have different dimension onto the suitable dimension of input value $\{η,ηx,ηy\}$ with `OceanLight.convertwave!`.

In this example, we will consider the case of flat surface elevation. Hence, $\{η,ηx,ηy\}$ is equal to the matrix of zeros. 

Once all the input variables are in place, `OceanLight.interface!` calculate the refraction of the light between two medium given surface elevation attribution and return the position, reflectance angle, and transmission ratio. 

```@example Center
OceanLight.interface!(xpb,ypb,zpb,θ,ϕ,fres,η,ηx,ηy,p)
```

`OceanLight` tracks the path of each photon travelling inside the water medium and store the irradiance value in the grid `ed`. 

Users need to specify these variables and corresponding dimension. 

```@example Center
ed = zeros(p.nx, p.ny, p.nz)
esol = zeros(p.num, p.nz)
randrng = MersenneTwister(1234)
area = zeros(4)
interi = zeros(Int64,4)
interj = zeros(Int64,4)
ix = div(p.nxη,2)+1
iy = div(p.nyη,2)+1
ϕps,θps = OceanLight.phasePetzold()
```

## Monte Carlo Simulation

`OceanLight` simulates the photon traveling inside the water medium, given its initial position $\{xpb,ypb,zpb\}$ and the direction it started with$\{θ,ϕ\}$. Once photons are inside the water, `OceanLight` will track its path, governed by its probability distribution and the attenuated coefficient input, and store the irradiance value in the grid `ed`. 

The `transfer!` function simulate a single photon path and store its landed position on the grid `ed`. Hence, to simulate multiple photons, users need to loop the `transfer!` function and giving the input of an individual photon's number `ip`. Thus, `OceanLight` could facilitate parallel computation. 

```@example Center
for ip = 1:p.nphoton
    OceanLight.transfer!(ed,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,xpb[ix,iy],
        ypb[ix,iy],zpb[ix,iy],area,interi,interj,randrng,η,ϕps,θps,p,1)
end
```

Lastly, once the field `ed` is obtained. The `applybc!` apply and ensure the boundary condition of the `ed`. 

```@example Center
OceanLight.applybc!(ed,p)
```

Once the simulation is over, the downwelling irradiance field $I(x,y,z)$ is stored in variable `ed`. To export the data, the build-in function `OceanLight.exported` export data and its statistical information in to `.h5` format, depended on the path given by user. 

The result downwelling irradiance field $I(x,y,z)$ is in 3 dimension tensor, where the first two dimensions represent horizontal field, and the last represents the vertical field or depth. The first few depth layer represents the air phase, hence $I(x,y,z)$ is zeros, depending on how many layer specified as air phases `ztop`. To visualize using `Plots` julia package, users can try the example code below. 

First, the irradiance field $I(x,y,z)$ is normalized by its maximum value $I_{0}$. To prevent undefined values `NaN` during logarithmic scaling, all zero-valued grid points are replaced with its minimum value.

```@example Center
max_val, max_loc = findmax(ed)
ed = ed./max_val
nonzero_vals = ed[ed .!= 0]
min_val = minimum(nonzero_vals)

for i in 1:Int(nxe+1)
    for j in 1:Int(nye+1)
        for k in ztop:nz
            if ed[i,j,k] == 0
                ed[i,j,k] = min_val
            end
        end
    end
end
```

```@example Center
using Plots
using Plots.Measures

# Choose slice indices
iy_c = 256   # y-index for vertical cross-section
iz_a = 40    # z-index for panel (a)
iz_b = 160   # z-index for panel (b)

# Define layout: 2 rows, 2 columns, but right column spans both rows
l = @layout [grid(2,1) c]

# Panel (a) : z = iz_a
p1 = heatmap(
    p.x .-10, p.y .-10, log.(ed[:,:,iz_a]),
    clim=(-20,-7), framestyle=:box, grid=false,
    c=cgrad(:viridis), legend=:none,
    xlabel="\$x(m)\$", ylabel="\$y(m)\$",
    title="(a) z = $(round(p.z[iz_a], digits=1)) m",
    xlim=[minimum(p.x).-10, maximum(p.x).-10],
    ylim=[minimum(p.y).-10, maximum(p.y).-10])
plot!(p1, [minimum(p.x).-10, maximum(p.x).-10], [p.y[iy_c]-10, p.y[iy_c]-10],
      color=:red, lw=2, ls=:dash, alpha=0.6)

# Panel (b) : z = iz_b
p2 = heatmap(
    p.x .-10, p.y .-10, log.(ed[:,:,iz_b]),
    clim=(-20,-7), framestyle=:box, grid=false,
    c=cgrad(:viridis), legend=:none,
    xlabel="\$x(m)\$", ylabel="\$y(m)\$",
    title="(b) z = $(round(p.z[iz_b], digits=1)) m",
    xlim=[minimum(p.x).-10, maximum(p.x).-10],
    ylim=[minimum(p.y).-10, maximum(p.y).-10])
plot!(p2, [minimum(p.x).-10, maximum(p.x).-10], [p.y[iy_c]-10, p.y[iy_c]-10],
      color=:red, lw=2, ls=:dash, alpha=0.6)

# Panel (c) : vertical cross-section at y = iy_c
p3 = heatmap(
    p.x .-10, reverse(p.z), reverse(transpose(log.(ed[:,iy_c,:]))),
    clim=(-20,-7), framestyle=:box, grid=false,
    c=cgrad(:viridis), ylim=(-(nz*dz-10),0),
    xlabel="\$x(m)\$", ylabel="\$z(m)\$",
    cbar_title="\$\\ln\\frac{I(x,y,z)}{I_{0}}\$",
    title="(c) y = $(round(p.y[iy_c]-10, digits=1)) m",
    xlim=[minimum(p.x).-10, maximum(p.x).-10])

# Add horizontal lines for z = iz_a, iz_b
plot!(p3, [minimum(p.x).-10, maximum(p.x).-10], [p.z[iz_a], p.z[iz_a]],
      color=:red, lw=1.5, ls=:dash, label="", alpha=0.6)
plot!(p3, [minimum(p.x).-10, maximum(p.x).-10], [p.z[iz_b], p.z[iz_b]],
      color=:red, lw=1.5, ls=:dash, label="", alpha=0.6)

# Combine
plot(p1, p2, p3, layout=l,
     size=(900,700),
     titleloc=:left, titlefont=font(8),
     left_margin=10mm, right_margin=10mm)

```

## Reference 

[^1]: Smith, R. C., & Baker, K. S. (1981). Optical properties of the clearest natural waters (200-800 nm). Applied optics, 20(2), 177–184. https://doi.org/10.1364/AO.20.000177 

[^2]: Kirk, J. T. O. (1981). Monte Carlo procedure for simulating the penetration of light into natural waters. In Technical paper - Commonwealth Scientific and Industrial Research Organization (Vol. 36).