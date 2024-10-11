# Air-Water Interaction

This section describes the function that can be used for the light interaction between the atmosphere and water. This section can be skipped entirely, if the simulation only considered the light within the water medium.

To setup the stage for light refraction, the data of water surface distribution is needed: water elevation `\eta`, the slope in x direction `\eta_{x} = \frac{\partial \eta}{\partial x}` and y direction `\eta_{y} = \frac{\partial \eta}{\partial y}`. OceanLight.jl offers random wave distribution or input in `.h5` format. 

## Setting random water surface distribution

```@docs
setwave!(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},rms::Float64,p::Param)
```

## Importing water surface distribution data

```@docs
readdata(datdir::String,fname::String,n::Tuple{<:Int64,<:Int64},pexy::Tuple{<:Float64,<:Float64})
```
The water surface distribution data may or may not match the input irradiance field. Hece, OceanLight.jl interpolates the value `\eta`, `\eta_x`, and `\eta_y`, to match the corresponding irradiance field. 

```@docs
convertwave!(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},
                      η0::Array{<:AbstractFloat,2},ηx0::Array{<:AbstractFloat,2},ηy0::Array{<:AbstractFloat,2},kbc=0::Int64)
```

## Refraction

```@docs
interface!(xpb::Array{<:Float64,2},ypb::Array{<:Float64,2},zpb::Array{<:Float64,2},
                    θ::Array{<:Float64,2},ϕ::Array{<:Float64,2},fres::Array{<:Float64,2},
                    η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},p::Param)
```