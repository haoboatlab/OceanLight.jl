# Light within Water

```@docs
transfer!(ed::Array{<:Float64,3},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,ip::Int64,
                   xpb::Float64,ypb::Float64,zpb::Float64,area::Vector{Float64},interi::Vector{Int64},
                   interj::Vector{Int64},randrng,η::Array{<:Float64,2},ph::Array{<:Float64,1},
                   θps::Array{<:Float64,1},p::Param,mode=0::Int64)
```

```@docs
transfer!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
                   edk::Array{<:Int64,1},count::Array{<:Int64,1},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,
                   ip::Int64,xpb::Float64,ypb::Float64,zpb::Float64,randrng,η::Array{<:Float64,2},
                   ph::Array{<:Float64,1},θps::Array{<:Float64,1},p::Param,mode=0::Int64)
```
