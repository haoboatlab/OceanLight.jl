# Air-Water Interaction


```@docs
setwave!(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},rms::Float64,p::Param)
```

## Importing water surface distribution data

```@docs
readdata(datdir::String,fname::String,n::Tuple{<:Int64,<:Int64},pexy::Tuple{<:Float64,<:Float64})
```

```@docs
convertwave!(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},
                      η0::Array{<:AbstractFloat,2},ηx0::Array{<:AbstractFloat,2},ηy0::Array{<:AbstractFloat,2},kbc=0::Int64)
```

## Refraction

```@docs
phasePetzold()
```

```@docs
interface!(xpb::Array{<:Float64,2},ypb::Array{<:Float64,2},zpb::Array{<:Float64,2},
                    θ::Array{<:Float64,2},ϕ::Array{<:Float64,2},fres::Array{<:Float64,2},
                    η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},p::Param)
```