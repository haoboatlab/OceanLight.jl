# Quick Start 
First, we install the OceanLight package.
```@example Center 
cd(mktempdir()) # hide

using Pkg 
Pkg.add("OceanLight")
using OceanLight 
```

##  Initial Condition 

All input variables required by OceanLight can be separated into 3 categories: Irradiance; setting up the Monte Carlo simulation grid field, 
Photon; number of photons and optical properties of the water, and Wave; structure of surface wave field. 

The detail on physical meaning of each input variables can be found here [Simulation.ModelSetup](@ref). 

```@example Center
# irradiance
nz = 200
dz = 1
nxe = 512
nye = 512
num = 31
ztop = 10
# photon
nphoton = 100000
kr = 10
nxp = 512
kbc = 0
b = 0.006
nyp = 512
a = 0.007
# wave
pey = 0.07853981633974483
nxeta = 512
nyeta = 512
pex = 0.07853981633974483
```
OceanLight accesses all input variable in the `Param` structure, throught the build-in function `readparams()` that access our input variable in `.yml` format.
The build-in function `writeparams()` will take the input of `dict` format and produce the `.yml` format. 

```@example Center
data=Dict("irradiance"=>Dict("nxe"=>nxe,"nye"=>nye,"nz"=>nz,"dz"=>dz,"ztop"=>ztop,"num"=>num),
            "wave"=>Dict("pex"=>pex,"pey"=>pey,"nxeta"=>nxeta,"nyeta"=>nyeta),
            "photon"=>Dict("nxp"=>nxp,"nyp"=>nyp,"nphoton"=>nphoton,"a"=>a,"b"=>b,"kr"=>kr,"kbc"=>kbc))

OceanLight.writeparams(data) 
p = OceanLight.readparams() 
```

## Air-Water Interaction

During the air-water interaction process, OceanLight simulates the photons transfer directly downward from the air side, interacts with the water surface, and transfer down into water medium. In addition to information inside structure `param`, users need to identify the surface elevation information: η; surface elevation, ηx; paritial derivation of η in x direction, and ηy; paritial derivation of η in y direction. 

The input value of η, ηx, and ηy need to have the same dimension with the incoming photons' grid. User can generate random surface elevation information with `OceanLight.setwave!`, or provided specific data. OceanLight can map the user's provided data of {η0,ηx0,ηy0} in different dimension onto the input value of {η,ηx,ηy} with `OceanLight.convertwave!`.
```@example Center
η = zeros(p.nxs,p.nys)
ηx = zeros(p.nxs,p.nys)
ηy = zeros(p.nxs,p.nys)
```

After photons' interaction with the surface, OceanLight requires the information of specific coordinate of photon in cartesian grid {xpb,ypb,zpb}, the direction in which photon will travel in polar coordinate {θ,ϕ}, and the fraction of light that transmit through the water {fres}: all in the dimension of incoming photon grid size. 
```@example Center
xpb = zeros(p.nxp,p.nyp)
ypb = zeros(p.nxp,p.nyp)
zpb = zeros(p.nxp,p.nyp)
θ = zeros(p.nxp,p.nyp)
ϕ = zeros(p.nxp,p.nyp)
fres = zeros(p.nxp,p.nyp)
```

With all the input and output variable specified, OceanLight calculate the interaction with `OceanLight.interface!`.
```@example Center
OceanLight.interface!(xpb,ypb,zpb,θ,ϕ,fres,η,ηx,ηy,p)
```

## Monte Carlo simulation

`OceanLight` simulates the photon traveling inside the water medium, given its initial position {xpb,ypb,zpb} and the direction it started with {θ,ϕ}. Once photons are inside the water, `OceanLight` will track its path, governed by its probability distribution and the attenuated coefficient input, and store the irradiance value in the grid `ed`. 

Users need to specify these variables and corresponding dimension. The detail on meaning of each input variables can be found [Simulation.WithinWater](@ref). 

```@example Center
using Random

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

The `transfer!` function simulate a single photon path and store its irradiance value on the grid `ed`. Hence, to simulate multiple photons, users need to loop the `transfer!` function and giving the input of an individual photon's number `ip`. Thus, `OceanLight` could facilitate parallel computation. 

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

## Export data

`OceanLight` exported the irradiance field `ed` and its statistics in `.h5` file. Detail on each specific input and the mode can be found [`exported`](@ref). 

```@example Center
OceanLight.exported(ed,η,p,"ed","3D",0)
```