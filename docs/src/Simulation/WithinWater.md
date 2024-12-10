# Light within water

In OceanLight.jl simulation, the scattering photons inside the water medium using the distribution based on Petzold(1972,) in which the value need to be called manually and stored in the variable `ph` and `θps`. 

```@docs
phasePetzold()
```
OceanLight.jl gives the options of constructing the irradiance solution field, depending on the initial condition and the problem one wish to solve. 

## Single point source  

If the problem is consisted of incoming photons on a single point (i.e. incoming photons only at the center of irradiance field,) the solution grid field is stored directly at `ed` solution field. 

```@docs
transfer!(ed::Array{<:Float64,3},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,ip::Int64,
                   xpb::Float64,ypb::Float64,zpb::Float64,area::Vector{Float64},interi::Vector{Int64},
                   interj::Vector{Int64},randrng,η::Array{<:Float64,2},ph::Array{<:Float64,1},
                   θps::Array{<:Float64,1},p::Param,mode=0::Int64)
```

## Multiple points source

If the problem is consisted of incoming photons on a multiple point (i.e. incoming photons everywhere on irradiance field,) OceanLight.jl will track each photon individually and stored the irradiance contribution on `ed1d` with its coordination at `edi`, `edj`, and `edk`, and later, combine into a single solution. 

```@docs
transfer!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
                   edk::Array{<:Int64,1},count::Array{<:Int64,1},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,
                   ip::Int64,xpb::Float64,ypb::Float64,zpb::Float64,randrng,η::Array{<:Float64,2},
                   ph::Array{<:Float64,1},θps::Array{<:Float64,1},p::Param,mode=0::Int64)
```
