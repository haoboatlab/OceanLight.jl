# Implementing MPI

```@setup  MPI_tutorial
cd(mktempdir()) 
using Pkg 
Pkg.add("OceanLight") 
Pkg.add("Plots") 
Pkg.add("MPI")
```

## Problem

```@example MPI_tutorial 
using OceanLight 
using Random
```
```@example MPI_tutorial 
# irradiance
nz = 200; dz = 1; nxe = 512; nye = 512; num = 31; ztop = 10                 
# photon
nphoton = 10000000; kr = 10; nxp = 512; kbc = 0; b = 0.0031; nyp = 512; a = 0.0196;                 
# wave
pey = 2*pi/20.0; nxeta = 512; nyeta = 512; pex = 2*pi/20.0;

data=Dict("irradiance"=>Dict("nxe"=>nxe,"nye"=>nye,"nz"=>nz,"dz"=>dz,"ztop"=>ztop,"num"=>num),
            "wave"=>Dict("pex"=>pex,"pey"=>pey,"nxeta"=>nxeta,"nyeta"=>nyeta),
            "photon"=>Dict("nxp"=>nxp,"nyp"=>nyp,"nphoton"=>nphoton,"a"=>a,"b"=>b,"kr"=>kr,"kbc"=>kbc))

OceanLight.writeparams(data)
```

## Initialize MPI

```@example MPI_tutorial 
using MPI
MPI.Init()
```

```@example MPI_tutorial 

```









```
mpirun -np 4 julia FILE_NAME.jl
```


