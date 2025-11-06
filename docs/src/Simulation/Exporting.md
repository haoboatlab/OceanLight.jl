# Exporting data

Once the HydrOptics.jl simulation completes, our solution grid field can be accessed from variable `ed`. 

## Periodic boundary condition

Before accessing the data, HydrOptics.jl offers a function in which applying Periodic Boundary Condition on the solution fied `ed`.

```@docs
applybc!(ed::Array{<:Float64,3},p::Param)
```

## Exporting data to HDF file

Our data can be exported to .h5 file, with 3 modes: `2D`, `3D`, and `full`. 
In `2D` mode, `.h5` file will store the statistic of our irradiance field solution; storing mean `μ`, variance `σ`, and cv `cv`, and the 2 dimenstion cross-section of our solution: `xz`, `yz`, and `xy`. 

If specified `3D` mode, in addition to the data from `2D` mode, HydrOptics.jl will stored the `ed` solution. 

If specified `full` mode, in addition to the data from `3D` mode, HydrOptics.jl will stored the physical coordination of the `ed` solution. 

```@docs
exported(ed::Array{<:Real,3},η::Array{<:Real,2},p::Param,
                  fname::String,mode="2D"::String,nk=0::Int64)
```

