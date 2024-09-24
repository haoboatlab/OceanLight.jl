# Model Setup

## Input Variable

```@docs
writeparams(data::Dict,fname="light.yml"::String)
```

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

## Simulation Parameters

```@docs
readparams(fname="light.yml"::String)
```

```@docs
Param    
```