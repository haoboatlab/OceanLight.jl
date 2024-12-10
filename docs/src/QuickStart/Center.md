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

The detail on physical meaning of each input variables can be found here (insert example block). 

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

OceanLight.writeparams() 
OceanLight.readparams() 
```

## Air-Water Interaction



```@example Center
η = zeros(p.nxs,p.nys)
ηx = zeros(p.nxs,p.nys)
ηy = zeros(p.nxs,p.nys)
```

```@example Center
xpb = zeros(p.nxp,p.nyp)
ypb = zeros(p.nxp,p.nyp)
zpb = zeros(p.nxp,p.nyp)
θ = zeros(p.nxp,p.nyp)
ϕ = zeros(p.nxp,p.nyp)
fres = zeros(zeros(p.nxp,p.nyp))
```

```@example Center
oceanlight.interface!(xpb,ypb,zpb,θ,ϕ,fres,η,ηx,ηy,p)
```

## Monte Carlo simulation

```@example Center
ed = zeros(p.nx, p.ny, pr.nz)
esol = zeros(p.num, p.nz)
randrng = MersenneTwister(1234)
area=zeros(4)
interi=zeros(Int64,4)
interj=zeros(Int64,4)
ix=div(p.nxη,2)+1
iy=div(p.nyη,2)+1

@time begin
    for ind=1:p.nphoton
        ip=allind[ind]
        OceanLight.transfer!(ed,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,xpb[ix,iy],
        ypb[ix,iy],zpb[ix,iy],area,interi,interj,randrng,η,ϕps,θps,parameter,1)
    end
end
```

```@example Center

```