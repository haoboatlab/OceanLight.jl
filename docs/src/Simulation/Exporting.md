# Exporting data

## Periodic Boundary Condition
```@docs
applybc!(ed::Array{<:Float64,3},p::Param)
```

## Exporting data to .h5 file

Our data can be exported to .h5 file, with 3 modes: `2D`, `3D`, and `full`. 
```@docs
exported(ed::Array{<:Real,3},η::Array{<:Real,2},p::Param,
                  fname::String,mode="2D"::String,nk=0::Int64)
```

