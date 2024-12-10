# Model setup

OceanLight.jl simulation reads all required variable from the struct `param`. The description on how to generate and setup the model is described below in this section.  

## Input variable

OceanLight.jl accepts and reads the input variables from .yml file. The structure on what variable to includes is shown below. 

```@docs
writeparams(data::Dict,fname="light.yml"::String)
```
This function only read the `dict` format. Hence, we first need to rearrange into `dict` format, before we can called the function `writeparams()`. 

```@example 
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

data=Dict("irradiance"=>Dict("nxe"=>nxe,"nye"=>nye,"nz"=>nz,"dz"=>dz,"ztop"=>ztop,"num"=>num),
            "wave"=>Dict("pex"=>pex,"pey"=>pey,"nxeta"=>nxeta,"nyeta"=>nyeta),
            "photon"=>Dict("nxp"=>nxp,"nyp"=>nyp,"nphoton"=>nphoton,"a"=>a,"b"=>b,"kr"=>kr,"kbc"=>kbc))
```

## Simulation parameters

OceanLight.jl reads all input variables in .yml format through the function `readparams()`, and store the values in structure `Param`. Beside our provided value, it will auto-generate some of the parameters that will be used in the simulation. The list and description of all the values can be accessed below. 

```@docs
readparams(fname="light.yml"::String)
```

```@docs
Param    
```