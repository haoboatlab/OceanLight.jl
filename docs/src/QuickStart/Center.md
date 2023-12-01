# Quick Start
```@example 1
tempdir() # hide 
```
## Initial Condition 

For the input initial condition, the function readparam will automatically read the input from light.yml file. 
Below, is the structure contained in the light.yml file.

```@example 1
using YAML

YAML.write_file("light.yml","
irradiance:
  nz: 200
  dz: 1
  nxe: 512
  nye: 512
  num: 31
  ztop: 10
photon:
  nphoton: 100000
  kr: 10
  nxp: 512
  kbc: 0
  b: 0.006
  nyp: 512
  a: 0.007
wave:
  pey: 0.07853981633974483
  nxeta: 512
  nyeta: 512
  pex: 0.07853981633974483
")
```
We, then, store all the input variable in the struct p, by using the function below. 
```@example 1
using Pkg # hide
Pkg.activate("https://github.com/tayT0T/LightMC.jl.git") # hide
import LightMC
p=readparams()
```
```@example 1
allind=1:p.nphoton
"η is the height(z axis) corresponding to each grid point on the surface, 2d array of grid number in x and y direction"
η=zeros(p.nxs,p.nys)
"ηx is the x coordination corresponding to each grid point on the surface, 2d array of grid number in x and y direction"
ηx=zeros(p.nxs,p.nys)
"ηy is the y coordination corresponding to each grid point on the surface, 2d array of grid number in x and y direction"
ηy=zeros(p.nxs,p.nys)
ϕps,θps=phasePetzold()
```

```@example 1
ed=zeros(p.nx,p.ny,p.nz)
esol=zeros(p.num,p.nz)
area=zeros(4)
interi=zeros(Int64,4)
interj=zeros(Int64,4)
xpb=zeros(p.nxp,p.nyp)
ypb=zeros(p.nxp,p.nyp)
zpb=zeros(p.nxp,p.nyp)
θ=zeros(p.nxp,p.nyp)
ϕ=zeros(p.nxp,p.nyp)
fres=zeros(p.nxp,p.nyp)
ix=div(p.nxη,2)+1
iy=div(p.nyη,2)+1
```


## Run the Monte Carlo Simulation
