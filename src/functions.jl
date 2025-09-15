include("SerialFFT.jl")
using .SerialFFT
using Test
using YAML 
using Random
using Statistics
using HDF5
using DocStringExtensions

"""
    struct Param    

All the Parameters and their attributes that will be used in the simulation
$(TYPEDFIELDS)
"""
struct Param    
    #Irradiance energy
    "nx = nxe+1 , nxe is the number of energy grid that we will calculate for the energy in x direction"
    nx::Int64
    "ny = nye+1 , nye is the number of energy grid that we will calculate for the energy in y direction"
    ny::Int64
    "number of energy grid that we will calculate for the energy in x direction"
    nxe::Int64 
    "number of energy grid that we will calculate for the energy in y direction"
    nye::Int64 
    "constant always set at 31 based on Measurement from Kirk, 1981 (The number of angle measurement in Kirk paper)"
    num::Int64
    "number of the layer beneath the sea surface (in z direction) that we will calculate for the irradiance energy"
    nz::Int64  
    "array of all the x corrdination"  
    x::Array{Float64,1}
    "array of all the y cooradination"
    y::Array{Float64,1}
    "array of all the z coordination"
    z::Array{Float64,1}
    #Domain
    "pex is being used in function pdfx or the function to calculate the partial derivative of the surface elevation in x direction"
    pex::Float64
    "pey is being used in function pdfx or the function to calculate the partial derivative of the surface elevation in y direction"
    pey::Float64
    "the distance between each grid point in x direction"
    dx::Float64
    "the distance between each grid point in y direction"
    dy::Float64
    "height difference between each layer that we calculate for the irradiance energy (nz)"
    dz::Float64
    "the maximum height above the water surface"
    ztop::Float64
    "minimum x value (always set to 0)"
    xmin::Float64
    "maximum x value equal to multiplication between total number of grid point (nxe) and the distance between each grid point (dx) in x direction"
    xmax::Float64
    "minimum y value (always set to 0)"
    ymin::Float64
    "maximum y value equal to multiplication between total number of grid point (nye) and the distance between each grid point (dy) in y direction"
    ymax::Float64
    "difference between the maximum and minimum value in x direction"
    xl::Float64
    "difference between the maximum and minimum value in y direction"
    yl::Float64
    #Surface elevation 
    "nxs=nxη+1"
    nxs::Int64 
    "nys=nyη+1"
    nys::Int64
    "same as nxeta value in the light.yml file: the number of wave grid point in x direction"
    nxη::Int64
    "same as nyeta value in the light.yml file: the number of wave grid point in y direction"
    nyη::Int64
    #Photon emission
    "number of grid in x direction that the photon will be emitted"
    nxp::Int64
    "number of grid in y direction that the photon will be emitted"
    nyp::Int64
    "number of photon emitted"
    nphoton::Int64
    "setting the mode of boundary condition: `kbc=1` (no interpolation) or periodic BC `kbc=0` (interpolation using FFT)"
    kbc::Int64
    "absortance coefficient"
    a::Float64
    "scattering coefficient"
    b::Float64
    "the multiple of sphere detector radius to dz (not being used anywhere in the code)"
    kr::Float64
    "the distance difference in x direction in which we will emit the photon (total length of physical grid xl devided by total number of photon grid nxp)"
    ddx::Float64
    "the distance difference in y direction in which we will emit the photon (total length of physical grid yl devided by total number of photon grid nyp)"
    ddy::Float64

    function Param(nxe::Integer,nye::Integer,num::Integer,nz::Integer,pex::Real,pey::Real,dz::Real,
                   ztop::Real,nxη::Integer,nyη::Integer,nxp::Integer,nyp::Integer,nphoton::Integer,
                   kbc::Integer,a::Real,b::Real,kr::Real)
        nx=nxe+1;         ny=nye+1
        "lenght difference between each grid point in x and y direction"
        dx=2*pi/pex/nxe;  dy=2*pi/pey/nye
        nxs=nxη+1;        nys=nyη+1
        xmin=0;           xmax=nxe*dx
        ymin=0;           ymax=nye*dy
        xl=xmax-xmin;     yl=ymax-ymin
        "an x and y grid where each point corresponding to the distance from (0,0)"
        x=[(i-1)*dx for i=1:nx]
        y=[(j-1)*dy for j=1:ny]
        "z: 1d array from the maximum value (in the atmosphere) to the deepest or the most negative z array"
        z=[-k*dz for k=1:nz] .+ ztop
        "the distance difference in which we will emit the photon (total length of the calculation grid/number of grid in which we want to emit photons)"
        ddx=xl/nxp;       ddy=yl/nyp

        new(nx,ny,nxe,nye,num,nz,x,y,z,pex,pey,dx,dy,dz,ztop,xmin,xmax,ymin,ymax,xl,yl,
            nxs,nys,nxη,nyη,nxp,nyp,nphoton,kbc,a,b,kr,ddx,ddy)
    end 
end

""" 
    readparams(fname="light.yml")

Read parameters from yml file `fname`. 

If `fname` is unspecified, use `light.yml` as the default file name.
"""
function readparams(fname="light.yml"::String)
    data = YAML.load_file(fname)

    key = "irradiance"
    nxe = data[key]["nxe"]
    nye = data[key]["nye"]
    nz = data[key]["nz"]
    dz = data[key]["dz"]
    ztop = data[key]["ztop"]
    num = data[key]["num"]

    key = "wave"
    pex = data[key]["pex"]
    pey = data[key]["pey"]
    nxη = data[key]["nxeta"]
    nyη = data[key]["nyeta"]
    
    key = "photon"
    nxp = data[key]["nxp"]
    nyp = data[key]["nyp"]
    nphoton = data[key]["nphoton"]
    kbc = data[key]["kbc"]
    a = data[key]["a"]
    b = data[key]["b"]
    kr = data[key]["kr"]
    
    @assert nxe > 1 "Invalid value of nxe"
    @assert nye > 1 "Invalid value of nye"
    @assert nz > 1 "Invalid value of nz"
    @assert dz > 0 "Invalid value of dz"
    @assert ztop > 0 "Invalid value of ztop"
    @assert num > 0 "Invalid value of num"
    @assert pex > 0 "Invalid value of pex"
    @assert pey > 0 "Invalid value of pey"
    @assert nxη > 1 "Invalid value of nxη" 
    @assert nyη > 1 "Invalid value of nyη"
    @assert nxp > 1 "Invalid value of nxp"
    @assert nyp > 1 "Invalid value of nyp"
    @assert nphoton >= 1 "Invalid value of nphoton"
    @assert typeof(kbc)<:Integer "Invalid value of kbc"
    @assert a > 0 "Invalid value of a"
    @assert b > 0 "Invalid value of b"
    @assert kr > 0 "Invalid value of kr"
    return Param(nxe,nye,num,nz,pex,pey,dz,ztop,nxη,nyη,nxp,nyp,nphoton,kbc,a,b,kr)
end

""" 
    writeparams(p::Param,fname="light.yml"::String)

Replace the data in the yml file `fname` to the new struct Param `P`.
    
If `fname` is unspecified, use `light.yml` as the default file name.
"""
function writeparams(p::Param,fname="light.yml"::String)    
    data=Dict("irradiance"=>Dict("nxe"=>p.nxe,"nye"=>p.nye,"nz"=>p.nz,"dz"=>p.dz,"ztop"=>p.ztop,"num"=>p.num),
              "wave"=>Dict("pex"=>p.pex,"pey"=>p.pey,"nxeta"=>p.nxη,"nyeta"=>p.nyη),
              "photon"=>Dict("nxp"=>p.nxp,"nyp"=>p.nyp,"nphoton"=>p.nphoton,"a"=>p.a,"b"=>p.b,"kr"=>p.kr))   
    writeparams(data,fname)
    return nothing
end

""" 
    writeparams(data::Dict,fname="light.yml"::String)

Replace the data in the yml file `fname` to the new dictionary `data`. 

If `fname` is unspecified, use `light.yml` as the default file name.
"""
function writeparams(data::Dict,fname="light.yml"::String)    
    YAML.write_file(fname,data)    
    return nothing
end

"""
    eta2rad!(ed::Array{<:Float64,3},esol::Array{<:Float64,2},η::Array{<:Float64,2},
                  ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},ph::Array{<:Float64,1},
                  θps::Array{<:Float64,1},p::Param,randrng)

Using the Monte Carlo to calculate the direction and angle of the photon that enter the water surface
"""
function eta2rad!(ed::Array{<:Float64,3},esol::Array{<:Float64,2},η::Array{<:Float64,2},
                  ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},ph::Array{<:Float64,1},
                  θps::Array{<:Float64,1},p::Param,randrng)
   # Calculate the photon that enter the wave surface, by using Monte Carlo Method, and also the periodic boundary conditon
    for ix=1:p.nxp, iy=1:p.nyp, ip=1:p.nphoton
        "Using the Monte Carlo Method"
        # photon emitted and enter water surface
        xpb,ypb,zpb,θ,ϕ,fres=interface(ix,iy,η,ηx,ηy,p)
        println("new photon emitted at i=$ix, j=$iy !")
        # photon travelling in medium
        transfer!(ed,esol,θ,ϕ,fres,ip,xpb,ypb,zpb,randrng,η,ph,θps,p,1)
    end
  
    if p.kbc == 0
        "Calculate the periodic boundary condition"
        for k=1:p.nz
            ed[1,1:p.nye,k]=ed[1,1:p.nye,k]+ed[p.nx,1:p.nye,k]
            ed[1:p.nxe,1,k]=ed[1:p.nxe,1,k]+ed[1:p.nxe,p.ny,k]
            ed[1,1,k]=ed[1,1,k]+ed[p.nx,p.ny,k]
        end
    end
end

function eta2rad!(ed::Array{<:Float64,3},esol::Array{<:Float64,2},η::Array{<:Float64,2},
                  ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},ph::Array{<:Float64,1},
                  θps::Array{<:Float64,1},ix::Int64,iy::Int64,ip::Int64,
                  p::Param,randrng)
    # Calculat only the photon that emitted and enter the wave surface, by using Monte Carlo Method
    xpb,ypb,zpb,θ,ϕ,fres=interface(ix,iy,η,ηx,ηy,p)
    # photon travelling in medium
    transfer!(ed,esol,θ,ϕ,fres,ip,xpb,ypb,zpb,randrng,
              η,ph,θps,p,1)
end

"""
    applybc!(ed::Array{<:Float64,3},p::Param)

Applying the periodic boundary condition on our irradiance output grid `ed` if `kbc` sets to 1 
"""
function applybc!(ed::Array{<:Float64,3},p::Param)
    if p.kbc == 0
        for k=1:p.nz
            ed[1,1:p.nye,k]=ed[1,1:p.nye,k]+ed[p.nx,1:p.nye,k]
            ed[1:p.nxe,1,k]=ed[1:p.nxe,1,k]+ed[1:p.nxe,p.ny,k]
            ed[1,1,k]=ed[1,1,k]+ed[p.nx,p.ny,k]
        end
    end
end

"""
    setwave!(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},rms::Float64,p::Param)

Giving the value to the existed wave surface distribution grid
# Arguments
- `η::Array{<:Float64,2}`: water surface elevation.
- `ηx::Array{<:Float64,2}`: partial derivative of water surface elevation in x direction.
- `ηy::Array{<:Float64,2}`: partial derivative of water surface elevation in y direction.
- `rms::Float64`: random number.
- `p::Param`: simulation parameters.
"""
function setwave!(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},rms::Float64,p::Param)
    η[:,:]=rand(p.nxη,p.nyη)    # /eta is the wave height function
    η[:,:]=η.-mean(η)
    σ=sqrt(mean(η.^2))
    η[:,:]=filtering(η.*rms/σ,0.3)
    ηx[:,:]=pdfx(η,p.pex)
    ηy[:,:]=pdfy(η,p.pey)
    return nothing
end

"""
    convertwave!(η, ηx, ηy, η0, ηx0, ηy0, kbc)   

Convert the surface wave geometry `η0`, `ηx0`, `ηy0` to  `η`, `ηx`, `ηy` with the same size as irradiance field. 
Can be used for nonperiodic BC`kbc=1` (no interpolation) or periodic BC `kbc=0` (interpolation using FFT)
"""
function convertwave!(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},
                      η0::Array{<:AbstractFloat,2},ηx0::Array{<:AbstractFloat,2},ηy0::Array{<:AbstractFloat,2},kbc=0::Int64)    
    if kbc==1
        "nonperiodic boundary condition (no interpolation)"
        η[1:size(η0,1),1:size(η0,2)]=η0
        ηx[1:size(η0,1),1:size(η0,2)]=ηx0
        ηy[1:size(η0,1),1:size(η0,2)]=ηy0        
    end
    if kbc==0
        "periodic boundary condition (with interpolation)"
        if size(η,1) > size(η0,1) && size(η,2) > size(η0,2)
            nh=size(η).-1
            η[1:nh[1],1:nh[2]]=padding(η0,nh)
            ηx[1:nh[1],1:nh[2]]=padding(ηx0,nh)
            ηy[1:nh[1],1:nh[2]]=padding(ηy0,nh)
        end
        η[1:nh[1],end]=η[1:nh[1],1]
        ηx[1:nh[1],end]=ηx[1:nh[1],1]
        ηy[1:nh[1],end]=ηy[1:nh[1],1]
        η[end,1:nh[2]]=η[1,1:nh[2]]
        ηx[end,1:nh[2]]=ηx[1,1:nh[2]]
        ηy[end,1:nh[2]]=ηy[1,1:nh[2]]
        η[end,end]=η[1,1]
        ηx[end,end]=ηx[1,1]
        ηy[end,end]=ηy[1,1]
    end
    return nothing
end

"""
    readdata(datdir::String,fname::String,n::Tuple{<:Int64,<:Int64},pexy::Tuple{<:Float64,<:Float64})

Reading the wave surface distribution data (surface elevation η and partial derivative of it in x and y ηx and ηy) from the .h5 file,
given the directory `datdir` and the file name `fname`
"""
function readdata(datdir::String,fname::String,n::Tuple{<:Int64,<:Int64},pexy::Tuple{<:Float64,<:Float64})
    fid=h5open(datdir*fname,"r")
    η=read("eta_hos")        
    ηx=read("ex_hos")
    ηy=read("ey_hos")
    close(fid)
    # ηx=pdfx(η,pex)
    # ηy=pdfy(η,pey)
    nxhos,nyhos=size(η)
    nx,ny=n
    pex,pey=pexy
    if nxhos>nx & nyhos>ny
        η=filtering(η,n)
        ηx=pdfx(η,pex)
    else
        η=padding(η,n)
        ηy=pdfy(η,pey)
    end        
    return η,ηx,ηy
end

"""
    interface!(xpb::Array{<:Float64,2},ypb::Array{<:Float64,2},zpb::Array{<:Float64,2},
                    θ::Array{<:Float64,2},ϕ::Array{<:Float64,2},fres::Array{<:Float64,2},
                    η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},p::Param)

Calculate the reflection and refraction of the photon or light ray that transmit from the atmosphere to the water.
# Arguments
- `xpb::Array{<:Float64,2}`: initial x coordination of the photon.
- `ypb::Array{<:Float64,2}`: initial y coordination of the photon.
- `zpb::Array{<:Float64,2}`: initial z coordination of the photon.
- `θ::Array{<:Float64,2}`: angle of the light ray relative to the z axis: polar angle.
- `ϕ::Array{<:Float64,2}`: angle of the light ray relative to the x axis: azimuthal angle.
- `fres::Array{<:Float64,2}`: fresnel coefficient or fractional transmission for unpolarized light
- `η::Array{<:Float64,2}`: water surface elevation.
- `ηx::Array{<:Float64,2}`: partial derivative of water surface elevation in x direction.
- `ηy::Array{<:Float64,2}`: partial derivative of water surface elevation in y direction.
- `p::Param`: simulation parameters.
"""
function interface!(xpb::Array{<:Float64,2},ypb::Array{<:Float64,2},zpb::Array{<:Float64,2},
                    θ::Array{<:Float64,2},ϕ::Array{<:Float64,2},fres::Array{<:Float64,2},
                    η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},p::Param)
    #loop through all the photon grid (nxp,nyp)
    
    for ix=1:p.nxp, iy=1:p.nyp
        xpb[ix,iy],ypb[ix,iy],zpb[ix,iy],θ[ix,iy],ϕ[ix,iy],fres[ix,iy]=interface(ix,iy,η,ηx,ηy,p)
    end
    return nothing
end


function interface(η::Array{<:Float64,2},ηx::Array{<:Float64,2},ηy::Array{<:Float64,2},p::Param)
    #setting the array dimension for the output: xpb, ypb, zpb, θ, and ϕ

    "the x-axis coordination of the Photon's position at the beginning"
    xpb=zeros(p.nxp,p.nyp)
    "the y-axis coordination of the Photon's position at the beginning"
    ypb=zeros(p.nxp,p.nyp)
    "the z-axis coordination of the Photon's position at the beginning"
    zpb=zeros(p.nxp,p.nyp)
    "the angle of the light ray relative to the z axis: polar angle"
    θ=zeros(p.nxp,p.nyp)
    "the angle of the light ray relative to the x axis: azimuthal angle"
    ϕ=zeros(p.nxp,p.nyp)
    "the fraction of light ray that transmit into the water"
    fres=zeros(p.nxp,p.nyp)
    for ix=1:p.nxp, iy=1:p.nyp
        xpb[ix,iy],ypb[ix,iy],zpb[ix,iy],θ[ix,iy],ϕ[ix,iy],fres[ix,iy]=interface(ix,iy,η,ηx,ηy,p)
    end
    return xpb,ypb,zpb,θ,ϕ,fres
end

function interface(ix::Int64,iy::Int64,η::Array{<:Float64,2},ηx::Array{<:Float64,2},
                   ηy::Array{<:Float64,2},p::Param)
    x=p.x
    y=p.y
    nx=p.nx
    ny=p.ny
    nz=p.nz
    ddx=p.ddx
    ddy=p.ddy
    # initialize parameters for refraction
    "na --refraction coefficient of air"
    na=1.
    "nw --refraction coefficient of water"
    nw=1.34
    #initial position of photon (xpb,ypb,zpb)      
    xpb=x[1]+ddx*(ix-1)
    ypb=y[1]+ddy*(iy-1)
    if xpb>p.xl || ypb>p.yl
        error("$ddx $ddy $ix $iy $xpb $ypb")
    end
    
    area,i,j=interpolation2d(xpb,ypb,p)
    η0=0
    ηx0=0
    ηy0=0
    for k=1:4
        η0+=area[k]*η[i[k],j[k]]/p.dx/p.dy
        ηx0+=area[k]*ηx[i[k],j[k]]/p.dx/p.dy
        ηy0+=area[k]*ηy[i[k],j[k]]/p.dx/p.dy
    end
    zpb=η0
    #apply refraction rule at air-water interface, using the Fresnel Reflectance Equation 
    "gama -- angle of reflection"
    gama=acos(1/(ηx0^2+ηy0^2+1)^0.5)               
    "gamap -- angle of trasmission"        
    gamap=asin(((ηx0^2+ηy0^2)/(ηx0^2+ηy0^2+1))^0.5*na/nw)  
    "θ is polar angle; angle between z axis and transmitted ray, when the photons is coming directly downward [0;0;-1]"
    θ=gama-gamap                
    if abs(gama) < 1e-12  
        "if the angle is near zero, all the energy will transmit into the water, without any reflection" 
        "all being trasmit, the fraction of trasmissiong is 1"
        fres=1.0                                      
        "azimuthal angle equals to 0"      
        ϕ=0.0
    else
        
        fres=((2*sin(gamap)*cos(gama)/sin(gama+gamap))^2
              +(2*sin(gamap)*cos(gama)/sin(gama+gamap)/cos(gama-gamap))^2)/2
        "temx is the normalizing term of the slope in the x direction, -1 <= ηx0 <=1; if the slope is negative, the light ray travel into the negative side, temx positive"
        temx=-ηx0/(ηx0^2+ηy0^2)^0.5
        "temy is the normalizing term of the slope in the y direction, -1 <= ηy0 <=1; if the slope is negative, the light ray travel into the negative side, temy positive "
        temy=-ηy0/(ηx0^2+ηy0^2)^0.5
        if abs(temy) <= 1e-12
            "If the surface plane does not have any slope in the y direction"
            if abs(temx-1) <= 1e-12
                "if normalized term of slope in x dir = -1, there is approx no slope in y direction: all the transmitted right ray going toward neg direction without any y coord"
                "azimuthal angle is pi; no y coordination"
                ϕ=pi
            end
            if abs(temx+1) <= 1e-12
                "if normalized term of slope in x dir = 1, there is approx no slope in y direction: all the transmitted right ray going toward pos direction without any y coord"
                "azimuthal angle is 0; no y coordination"
                ϕ=0
            end
        else 
            if temy > 0
                "if normalized term of slope in y dir is negative, y coord of scattered ray is negative, light in third and fourth quadrant"
                "if slope is positive, light ray is in positive x coord, temx is negative, arccos(temx) is higher than 90 deg, and azimuthal angle land in fourth quadrant"
                ϕ=pi+acos(temx)
            else 
                if temy < 0
                    "if normalized term of slope in y dir is positive, y coord of scattered ray is positive, light in first and second quadrant"
                    "if slope is positive, light ray is in positive x coord, temx is negative, arccos(temx) is higher than 90 deg, and azimuthal angle land in first quadrant"
                    ϕ=pi-acos(temx)
                end
            end
        end
    end      
    if (@isdefined ϕ)==false
        println("$xpb $ypb $zpb $η0 $ηx0 $ηy0")
        println("$gama $gamap $θ")
    end
    return xpb,ypb,zpb,θ,ϕ,fres
end   

"""
    transfer!(ed::Array{<:Float64,3},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,ip::Int64,
                   xpb::Float64,ypb::Float64,zpb::Float64,area::Vector{Float64},interi::Vector{Int64},
                   interj::Vector{Int64},randrng,η::Array{<:Float64,2},ph::Array{<:Float64,1},
                   θps::Array{<:Float64,1},p::Param,mode=0::Int64)

Doing the Monte Carlo Simulation.
# Arguments
- `ed::Array{<:Float64,3}`: Irradiance solution grid 
- `esol::Array{<:Float64,2}`: Irradiance solution grid for `solar` mode (under deverlopment)
- `θ::Float64`: angle of the light ray relative to the z axis: polar angle.
- `ϕ::Float64`: angle of the light ray relative to the x axis: azimuthal angle.
- `fres::Float64`: fresnel coefficient or fractional transmission for unpolarized light
- `ip::Int64`: current photon's number being simulated (ie. ip ∈ {1, 2,..., nphoton}) 
- `xpb::Float64`: initial x coordination of the photon.
- `ypb::Float64`: initial y coordination of the photon.
- `zpb::Float64`: initial z coordination of the photon.
- `area::Vector{Float64}`: 4 values of the area inside a single grid where a photon lands corresponding to: 4 corners of the square grid.
- `interi::Vector{Int64}`: x coordination (grid number) of a single grid where a photon lands: from bottom left, bottom right, upper left, and upper right.
- `interj::Vector{Int64}`: y coordination (grid number) of a single grid where a photon lands: from bottom left, bottom right, upper left, and upper right.
- `randrng`: PRNGs (pseudorandom number generators) exported by the Random package.
- `η::Array{<:Float64,2}`: water surface elevation.
- `ph::Array{<:Float64,1}`: cumulation distribution of scattering angle (obtained from `phasePetzold()`)
- `θps::Array{<:Float64,1}`: angle between new trajectory and the direction of the photon before scattering corresponding to each `ϕps` (obtained from `phasePetzold()`)
- `p::Param`: simulation parameters.
- `mode::Int64`: mode of different irradinace calculation (under deverlopment)
"""
function transfer!(ed::Array{<:Float64,3},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,ip::Int64,
                   xpb::Float64,ypb::Float64,zpb::Float64,area::Vector{Float64},interi::Vector{Int64},
                   interj::Vector{Int64},randrng,η::Array{<:Float64,2},ph::Array{<:Float64,1},
                   θps::Array{<:Float64,1},p::Param,mode=0::Int64)
    # initialize parameters
    nums=size(ph,1)
    x=p.x
    y=p.y
    z=p.z
    nx=p.nx
    ny=p.ny
    nz=p.nz
    dz=p.dz
    xl=p.xl
    yl=p.yl
    xmin=p.xmin
    xmax=p.xmax
    ymin=p.ymin
    ymax=p.ymax
    "beam attenuation coefficient: c = a + b"
    c=p.a+p.b
    "the variable determine the end of the photon; if idie==1, the photon is being absorb, and the montecarlo method will stop running"
    idie=0
    "the variable determine the scattering of the photon; if isca==1, the photon is scattered and will be continue in the loop"
    isca=0
    θs=0.0
    ϕs=0.0
    
    while true
        "Doing the Monte Carlo procedure, finding the energy level in each interaction, and looping until all the energy has been absorbed"
        "lint:the length photon travel until next interaction, dzl:projection of lint to vertical direction, inter: energy level grid that photon travel to"
        lint,dzl,inter=interlength(dz,θ,c,randrng)      
#       println((lint,dzl,inter))
        if inter == 0
            "if the photon end up at the surface level"
            "x coordination of where photon end up to: initial position of photon in x coord plus length photon travel in x direction"
            xpe=xpb+lint*sin(θ)*cos(ϕ)
            "y coordination of where photon end up to: initial position of photon in y coord plus length photon travel in y direction"
            ype=ypb+lint*sin(θ)*sin(ϕ)
            "z coordination of where photon end up to: initial position of photon in z coord minus (depth is in negative z direction) length photon travel in z direction"
            zpe=zpb-lint*cos(θ)  
            "number of the energy layer that the photons travel, from the top ztop"        
            kk=floor(Int,(z[1]-zpb)/dz)
            if dzl > 0
                "if the photon travel to some depth"
                "isign determine which direction photons travle, if isign ==1 photons travel downward, isign == -1, photons reflect and travel upward"
                isign=1
                kk=kk+1
            else 
                if dzl < 0
                    "if photon reflect back"
                    isign=-1
                end
            end
            "we normally work with mode 1 skip the solar_energy function"        
            if mode==2 || mode==0
                solar_energy!(esol,kk,θ,xpb,ypb,zpb,xpe,ype,zpe,isign,p)
            end
            if isign*(z[kk]-zpe) > 0
                "if the photon travle upward to some energy level on the surface calculate the energy at the ending point"
                if mode==1 || mode==0
                    energy!(ed,kk,fres,θ,area,interi,interj,xpb,ypb,zpb,xpe,ype,zpe,p)
                end
            end            
            if xpe < xmin||xpe>xmax||ype<ymin||ype>ymax
                "if the photon ends up traveling outside the grid we are interested in"
                if p.kbc==1
                    idie=1
                    break
                else 
                    if p.kbc==0
                        xpe,ype=regulatexy(xpe,ype,p)
                    end
                end
            end
            if isign == -1
                idie,isca=surface(idie,isca,xpe,ype,zpe,area,interi,interj,η,ip,p)
                if idie == 1
                    break
                end
            end 
            "at the surface calculate whether the photon reflect back to the water or refraction to the air"           
            θs,ϕs,isca,idie=interaction(p.a,p.b,ph,θps,randrng)             
            if idie == 1            
#                println("photon $ip dies at level $kk")
                break            
            else 
                if isca==1
                    θ,ϕ,θs,ϕs=scatter(θ,ϕ,θs,ϕs)
                    continue           
                else
                    println("Interaction error")
                    println("idie=$idie")
                    println("isca=$isca")
                end
            end
            xpb=xpe
            ypb=ype
            zpb=zpe
        else
            "if the photon travle downward (from the initial position at the water surface)"
            kk=floor(Int,(z[1]-zpb)/dz)
            if dzl>0
                isign=1
                kk=kk+1
            else
                isign=-1
            end        
            dl=abs(dz/cos(θ))
            # dl is a scalar, dl>0        
            for iint=1:inter 
                "looping calculate the energy level at every single level that photon piercing through"           
                if iint==inter
                    dl=dl-mod(lint,dl)
                end            
                xpe=xpb+dl*sin(θ)*cos(ϕ)
                ype=ypb+dl*sin(θ)*sin(ϕ)
                zpe=zpb-dl*cos(θ)            
                if mode==1 || mode==0
                    energy!(ed,kk,fres,θ,area,interi,interj,xpb,ypb,zpb,xpe,ype,zpe,p)
                end
                if mode==2 || mode==0
                    solar_energy!(esol,kk,θ,xpb,ypb,zpb,xpe,ype,zpe,isign,p)
                end
                if xpe < xmin||xpe>xmax||ype<ymin||ype>ymax
                    if p.kbc==1
                        idie=1
                        break
                    else 
                        if p.kbc==0
                            xpe,ype=regulatexy(xpe,ype,p)
                        end
                    end
                end
                if isign==-1
                    idie,isca=surface(idie,isca,xpe,ype,zpe,area,interi,interj,η,ip,p)
                    if idie == 1
                        break
                    end
                end
                kk=kk+isign
                xpb=xpe
                ypb=ype
                zpb=zpe
                if zpb < z[nz]
                    idie=1
#                    println("photon $ip reaches bottom and dies !")
                    break
                end
            end
            if idie==1
                break
            end
            "the first interaction between the photon and particle in the water"
            θs,ϕs,isca,idie=interaction(p.a,p.b,ph,θps,randrng)
            if idie==1         
#                println("photon $ip dies at level $kk")
                break         
            else 
                if isca == 1
                    θ,ϕ,θs,ϕs=scatter(θ,ϕ,θs,ϕs)
                    continue                          
                else
                    println("idie=$idie")
                    println("isca=$isca")
                    error("Interaction error")
                end
            end
        end
    end
    return nothing
end 

"""
    transfer!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
                   edk::Array{<:Int64,1},count::Array{<:Int64,1},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,
                   ip::Int64,xpb::Float64,ypb::Float64,zpb::Float64,randrng,η::Array{<:Float64,2},
                   ph::Array{<:Float64,1},θps::Array{<:Float64,1},p::Param,mode=0::Int64)

Doing the Monte Carlo Simulation.
# Arguments
- `ed1d::Array{<:Float64,1}`: fraction of irradiance that will be assigned to 4 corners of a grid where a photon lands
- `edi::Array{<:Int64,1}`: x coordination (grid number) of a single grid where a photon lands: from bottom left, bottom right, upper left, and upper right.
- `edj::Array{<:Int64,1}`: y coordination (grid number) of a single grid where a photon lands: from bottom left, bottom right, upper left, and upper right.
- `edk::Array{<:Int64,1}`: number of the energy layer that the photons travel, from the top ztop in 1 by 4 array
- `count::Array{<:Int64,1}`: (parrallel computing) dummy integer span from 1 to 4 to keep track of the size of the ed1d, edi, edj, edk. 
- `esol::Array{<:Float64,2}`:  Irradiance solution grid for `solar` mode (under deverlopment)
- `θ::Float64`: angle of the light ray relative to the z axis: polar angle.
- `ϕ::Float64`: angle of the light ray relative to the x axis: azimuthal angle.
- `fres::Float64`: fresnel coefficient or fractional transmission for unpolarized light.
- `ip::Int64`: current photon's number being simulated (ie. ip ∈ {1, 2,..., nphoton}) 
- `xpb::Float64`: initial x coordination of the photon.
- `ypb::Float64`: initial y coordination of the photon.
- `zpb::Float64`: initial z coordination of the photon.
- `randrng`: PRNGs (pseudorandom number generators) exported by the Random package
- `η::Array{<:Float64,2}`: water surface elevation.
- `ph::Array{<:Float64,1}`: cumulation distribution of scattering angle (obtained from `phasePetzold()`)
- `θps::Array{<:Float64,1}`: angle between new trajectory and the direction of the photon before scattering corresponding to each `ϕps` (obtained from `phasePetzold()`)
- `p::Param`: simulation parameters.
- `mode::Int64`: mode of different irradinace calculation (under deverlopment)
"""
function transfer!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
                   edk::Array{<:Int64,1},count::Array{<:Int64,1},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,
                   ip::Int64,xpb::Float64,ypb::Float64,zpb::Float64,randrng,η::Array{<:Float64,2},
                   ph::Array{<:Float64,1},θps::Array{<:Float64,1},p::Param,mode=0::Int64)
    # initialize parameters
    nums=size(ph,1)
    x=p.x
    y=p.y
    z=p.z
    nx=p.nx
    ny=p.ny
    nz=p.nz
    dz=p.dz
    xl=p.xl
    yl=p.yl
    xmin=p.xmin
    xmax=p.xmax
    ymin=p.ymin
    ymax=p.ymax
    c=p.a+p.b
    idie=0
    isca=0
    θs=0.0
    ϕs=0.0
    area=zeros(4)
    interi=zeros(Int64,4)
    interj=zeros(Int64,4)
    # Monte Carlo procedure    
    while true
        lint,dzl,inter=interlength(dz,θ,c,randrng)      # lint is the length before the photon and medium interact
#        println((lint,dzl,inter))
        if inter == 0
            xpe=xpb+lint*sin(θ)*cos(ϕ)
            ype=ypb+lint*sin(θ)*sin(ϕ)
            zpe=zpb-lint*cos(θ)          
            kk=floor(Int,(z[1]-zpb)/dz)
            if dzl > 0
                isign=1
                kk=kk+1
            else 
                if dzl < 0
                    isign=-1
                end
            end        
            if mode==2 || mode==0
                solar_energy!(esol,kk,θ,xpb,ypb,zpb,xpe,ype,zpe,isign,p)
            end
            if isign*(z[kk]-zpe) > 0
                if mode==1 || mode==0
                    energy!(ed1d,edi,edj,edk,count,kk,fres,θ,xpb,ypb,zpb,xpe,ype,zpe,p)
                end
            end            
            if xpe < xmin||xpe>xmax||ype<ymin||ype>ymax
                if p.kbc==1
                    idie=1
                    break
                else 
                    if p.kbc==0
                        xpe,ype=regulatexy(xpe,ype,p)
                    end
                end
            end
            if isign == -1
                idie,isca=surface(idie,isca,xpe,ype,zpe,area,interi,interj,η,ip,p)
                if idie == 1
                    break
                end
            end            
            θs,ϕs,isca,idie=interaction(p.a,p.b,ph,θps,randrng)             
            if idie == 1            
#                println("photon $ip dies at level $kk")
                break            
            else 
                if isca==1
                    θ,ϕ,θs,ϕs=scatter(θ,ϕ,θs,ϕs)
                    continue           
                else
                    println("Interaction error")
                    println("idie=$idie")
                    println("isca=$isca")
                end
            end
            xpb=xpe
            ypb=ype
            zpb=zpe
        else
            kk=floor(Int,(z[1]-zpb)/dz)
            if dzl>0
                isign=1
                kk=kk+1
            else
                isign=-1
            end        
            dl=abs(dz/cos(θ))
            # dl is a scalar, dl>0        
            for iint=1:inter            
                if iint==inter
                    dl=dl-mod(lint,dl)
                end            
                xpe=xpb+dl*sin(θ)*cos(ϕ)
                ype=ypb+dl*sin(θ)*sin(ϕ)
                zpe=zpb-dl*cos(θ)            
                if mode==1 || mode==0
                    energy!(ed1d,edi,edj,edk,count,kk,fres,θ,xpb,ypb,zpb,xpe,ype,zpe,p)
                end
                if mode==2 || mode==0
                    solar_energy!(esol,kk,θ,xpb,ypb,zpb,xpe,ype,zpe,isign,p)
                end
                if xpe < xmin||xpe>xmax||ype<ymin||ype>ymax
                    if p.kbc==1
                        idie=1
                        break
                    else 
                        if p.kbc==0
                            xpe,ype=regulatexy(xpe,ype,p)
                        end
                    end
                end
                if isign==-1
                    idie,isca=surface(idie,isca,xpe,ype,zpe,area,interi,interj,η,ip,p)
                    if idie == 1
                        break
                    end
                end
                kk=kk+isign
                xpb=xpe
                ypb=ype
                zpb=zpe
                if zpb < z[nz]
                    idie=1
#                    println("photon $ip reaches bottom and dies !")
                    break
                end
            end
            if idie==1
                break
            end
            θs,ϕs,isca,idie=interaction(p.a,p.b,ph,θps,randrng)
            if idie==1         
#                println("photon $ip dies at level $kk")
                break         
            else 
                if isca == 1
                    θ,ϕ,θs,ϕs=scatter(θ,ϕ,θs,ϕs)
                    continue                          
                else
                    println("idie=$idie")
                    println("isca=$isca")
                    error("Interaction error")
                end
            end
        end
    end
    return nothing
end 

# function transfer!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
#                    edk::Array{<:Int64,1},esol::Array{<:Float64,2},θ::Float64,ϕ::Float64,fres::Float64,
#                    ip::Int64,xpb::Float64,ypb::Float64,zpb::Float64,randrng,η::Array{<:Float64,2},
#                    ph::Array{<:Float64,1},θps::Array{<:Float64,1},p::Param,mode=0::Int64)
#     @warn "Deprecated low performance function!"

#     # initialize parameters
#     nums=size(ph,1)
#     x=p.x
#     y=p.y
#     z=p.z
#     nx=p.nx
#     ny=p.ny
#     nz=p.nz
#     dz=p.dz
#     xl=p.xl
#     yl=p.yl
#     xmin=p.xmin
#     xmax=p.xmax
#     ymin=p.ymin
#     ymax=p.ymax
#     c=p.a+p.b
#     idie=0
#     isca=0
#     θs=0.0
#     ϕs=0.0
#     # Monte Carlo procedure    
#     while true
#         lint,dzl,inter=interlength(dz,θ,c,randrng)      
# #        println((lint,dzl,inter))
#         if inter == 0
#             xpe=xpb+lint*sin(θ)*cos(ϕ)
#             ype=ypb+lint*sin(θ)*sin(ϕ)
#             zpe=zpb-lint*cos(θ)          
#             kk=floor(Int,(z[1]-zpb)/dz)
#             if dzl > 0
#                 isign=1
#                 kk=kk+1
#             else 
#                 if dzl < 0
#                     isign=-1
#                 end
#             end        
#             if mode==2 || mode==0
#                 solar_energy!(esol,kk,θ,xpb,ypb,zpb,xpe,ype,zpe,isign,p)
#             end
#             if isign*(z[kk]-zpe) > 0
#                 if mode==1 || mode==0
#                     energy!(ed1d,edi,edj,edk,kk,fres,θ,xpb,ypb,zpb,xpe,ype,zpe,p)
#                 end
#             end            
#             if xpe < xmin||xpe>xmax||ype<ymin||ype>ymax
#                 if p.kbc==1
#                     idie=1
#                     break
#                 else 
#                     if p.kbc==0
#                         xpe,ype=regulatexy(xpe,ype,p)
#                     end
#                 end
#             end
#             if isign == -1
#                 idie,isca=surface(idie,isca,xpe,ype,zpe,area,interi,interj,η,ip,p)
#                 if idie == 1
#                     break
#                 end
#             end            
#             θs,ϕs,isca,idie=interaction(p.a,p.b,ph,θps,randrng)             
#             if idie == 1            
# #                println("photon $ip dies at level $kk")
#                 break            
#             else 
#                 if isca==1
#                     θ,ϕ,θs,ϕs=scatter(θ,ϕ,θs,ϕs)
#                     continue           
#                 else
#                     println("Interaction error")
#                     println("idie=$idie")
#                     println("isca=$isca")
#                 end
#             end
#             xpb=xpe
#             ypb=ype
#             zpb=zpe
#         else
#             kk=floor(Int,(z[1]-zpb)/dz)
#             if dzl>0
#                 isign=1
#                 kk=kk+1
#             else
#                 isign=-1
#             end        
#             dl=abs(dz/cos(θ))
#             # dl is a scalar, dl>0        
#             for iint=1:inter            
#                 if iint==inter
#                     dl=dl-mod(lint,dl)
#                 end            
#                 xpe=xpb+dl*sin(θ)*cos(ϕ)
#                 ype=ypb+dl*sin(θ)*sin(ϕ)
#                 zpe=zpb-dl*cos(θ)            
#                 if mode==1 || mode==0
#                     energy!(ed1d,edi,edj,edk,kk,fres,θ,xpb,ypb,zpb,xpe,ype,zpe,p)
#                 end
#                 if mode==2 || mode==0
#                     solar_energy!(esol,kk,θ,xpb,ypb,zpb,xpe,ype,zpe,isign,p)
#                 end
#                 if xpe < xmin||xpe>xmax||ype<ymin||ype>ymax
#                     if p.kbc==1
#                         idie=1
#                         break
#                     else 
#                         if p.kbc==0
#                             xpe,ype=regulatexy(xpe,ype,p)
#                         end
#                     end
#                 end
#                 if isign==-1
#                     idie,isca=surface(idie,isca,xpe,ype,zpe,area,interi,interj,η,ip,p)
#                     if idie == 1
#                         break
#                     end
#                 end
#                 kk=kk+isign
#                 xpb=xpe
#                 ypb=ype
#                 zpb=zpe
#                 if zpb < z[nz]
#                     idie=1
#                     println("photon $ip reaches bottom and dies !")
#                     break
#                 end
#             end
#             if idie==1
#                 break
#             end
#             θs,ϕs,isca,idie=interaction(p.a,p.b,ph,θps,randrng)
#             if idie==1         
# #                println("photon $ip dies at level $kk")
#                 break         
#             else 
#                 if isca == 1
#                     θ,ϕ,θs,ϕs=scatter(θ,ϕ,θs,ϕs)
#                     continue                          
#                 else
#                     println("idie=$idie")
#                     println("isca=$isca")
#                     error("Interaction error")
#                 end
#             end
#         end
#     end
#     return nothing
# end 

"""
    regulatexy(xpe::Float64,ype::Float64,p::Param)

Applying the periodic boundary condition if the landed photons coordination (xpe and ype) is higher or lower 
than its boundary (xmax,ymax,xmin,xmax)
"""
function regulatexy(xpe::Float64,ype::Float64,p::Param)
    while xpe < p.xmin
        xpe=xpe+p.xl
    end
    while xpe > p.xmax
        xpe=xpe-p.xl
    end
    while ype < p.ymin
        ype=ype+p.yl
    end
    while ype > p.ymax
        ype=ype-p.yl
    end
    return xpe,ype
end

function regulateij(i::Array{<:Int64,1},j::Array{<:Int64,1},nx::Int64,ny::Int64)
    if i[1]==0
        i[1] += nx-1
        i[2] += nx-1
        i[3] += nx-1
        i[4] += nx-1
    end
    if i[4]==nx+1
        i[1] -= nx-1
        i[2] -= nx-1
        i[3] -= nx-1
        i[4] -= nx-1
    end
    if j[1]==0
        j[1] += ny-1
        j[2] += ny-1
        j[3] += ny-1
        j[4] += ny-1
    end
    if j[4]==ny+1
        j[1] -= ny-1
        j[2] -= ny-1
        j[3] -= ny-1
        j[4] -= ny-1
    end
    return i,j
end

"""
    regulateij!(i, j, nx, ny)

if the position exceed the boundary, that position equal to another end. 
Periodic Boundary condition (?), one end is equal to another end. 
"""
function regulateij!(i::Array{<:Int64,1},j::Array{<:Int64,1},nx::Int64,ny::Int64)
    if i[1]==0
        i[1] += nx-1
        i[2] += nx-1
        i[3] += nx-1
        i[4] += nx-1
    end
    if i[4]==nx+1
        i[1] -= nx-1
        i[2] -= nx-1
        i[3] -= nx-1
        i[4] -= nx-1
    end
    if j[1]==0
        j[1] += ny-1
        j[2] += ny-1
        j[3] += ny-1
        j[4] += ny-1
    end
    if j[4]==ny+1
        j[1] -= ny-1
        j[2] -= ny-1
        j[3] -= ny-1
        j[4] -= ny-1
    end
    return nothing
end
"""
    interlength(dz::Float64,θ::Float64,c::Float64,randrng)  

Calculate for the projection length of the photon after the scattering to another interaction, by using monte carlo method.
When c is the beam attenuate coefficient; c = absorbtance coeff(a) + scatterance coeff(b)
"""
function interlength(dz::Float64,θ::Float64,c::Float64,randrng)            
    "random the number to be used in the Monte Carlo Simulation"
    rad=rand(randrng)
    "lint is the length for interaction between photon and medium to happen"
    lint=-log(rad)/c      
    "dzl is the projection of lint to vertical direction"
    dzl=lint*cos(θ)
    "Since we need to calculate the energy in each specific layer nz, we round up the dzl, so that it correspond to our set up grid"
    inter=floor(Int,abs(dzl)/dz)    
    return lint,dzl,inter
end

"""
    energy!(ed::Array{<:Float64,3},kk::Int64,fres::Float64,θ::Float64,
                 area::Vector{Float64},interi::Vector{Int64},interj::Vector{Int64},
                 xpb::Float64,ypb::Float64,zpb::Float64,xpe::Float64,ype::Float64,zpe::Float64,p::Param)

Update the 3D irradiance field `ed'.   
# Arguments
- `ed::Array{<:Float64,3}`: Irradiance solution grid 
- `kk::Int64`: number of the energy layer that the photons travel, from the top ztop
- `fres::Float64`: Fresnel coefficient.
- `θ::Float64`: photon propagation direction.
- `area::Vector{Float64}`: 4 values of the area inside a single grid where a photon lands corresponding to: 4 corners of the square grid.
- `interi::Vector{Int64}`: x coordination (grid number) of a single grid where a photon lands: from bottom left, bottom right, upper left, and upper right.
- `interj::Vector{Int64}`: y coordination (grid number) of a single grid where a photon lands: from bottom left, bottom right, upper left, and upper right.
- `xpb::Float64`: initial x coordination of the photon.
- `ypb::Float64`: initial y coordination of the photon.
- `zpb::Float64`: initial z coordination of the photon.
- `xpe::Float64`: x coordination of where photon end up to
- `ype::Float64`: y coordination of where photon end up to
- `zpe::Float64`: z coordination of where photon end up to
- `p::Param`: simulation parameters                
"""
function energy!(ed::Array{<:Float64,3},kk::Int64,fres::Float64,θ::Float64,
                 area::Vector{Float64},interi::Vector{Int64},interj::Vector{Int64},
                 xpb::Float64,ypb::Float64,zpb::Float64,xpe::Float64,ype::Float64,zpe::Float64,p::Param)
    # check if downwelling
    if θ < pi/2
        # obtain crossing point on certain depth
        fract=(p.z[kk]-zpb)/(zpe-zpb)
        xpk=xpb+(xpe-xpb)*fract
        ypk=ypb+(ype-ypb)*fract    
        if xpk<p.xmin || xpk>=p.xmax || ypk<p.ymin || ypk>=p.ymax
            xpk,ypk=regulatexy(xpk,ypk,p)
        end
        
        # distribute ed to grid points
        interpolation2d!(area,interi,interj,xpk,ypk,p)
        # distribute energy to grid point on level kk
        for k=1:4
            ed[interi[k],interj[k],kk]+=(area[k]/p.dx/p.dy)*fres*cos(θ)
        end
    end
    return nothing
end

# """
# function energy!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
#                  edk::Array{<:Int64,1},kk::Int64,fres::Float64,θ::Float64,xpb::Float64,
#                  ypb::Float64,zpb::Float64,xpe::Float64,ype::Float64,zpe::Float64,p::Param)

# Input `kk' the vertical level, `fres' the Fresnel coefficient, `θ' the photon propagation direction, 
# (`xpb',`ypb',`zpb') the coordinates of the photon, (`xpe',`ype',`zpe') the coordinates of the irradiance field,
# and `p' the parameters.

# Output the irradiance strength `ed1d' and the information of the index `edi', `edj', and `edk', all as 1D arrays. 

# """
# function energy!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
#                  edk::Array{<:Int64,1},kk::Int64,fres::Float64,θ::Float64,xpb::Float64,
#                  ypb::Float64,zpb::Float64,xpe::Float64,ype::Float64,zpe::Float64,p::Param)
#     @warn "Deprecated low performance function!"
#     # check if downwelling
#     if θ < pi/2
#         # obtain crossing point on certain depth
#         fract=(p.z[kk]-zpb)/(zpe-zpb)
#         xpk=xpb+(xpe-xpb)*fract
#         ypk=ypb+(ype-ypb)*fract    
#         if xpk<p.xmin || xpk>=p.xmax || ypk<p.ymin || ypk>=p.ymax
#             xpk,ypk=regulatexy(xpk,ypk,p)
#         end
        
#         # distribute ed to grid points
#         area,i,j=interpolation2d(xpk,ypk,p)
#         # distribute energy to grid point on level kk
#         for k=1:4
#             push!(ed1d,(area[k]/p.dx/p.dy)*fres*cos(θ))
#             push!(edi,i[k])
#             push!(edj,j[k])
#             push!(edk,kk)
#             #            ed[i[k],j[k],kk]+=(area[k]/p.dx/p.dy)*fres*cos(θ)
#         end
#     end
#     return nothing
# end

"""
    energy!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
                 edk::Array{<:Int64,1},count::Array{<:Int64,1},kk::Int64,fres::Float64,θ::Float64,
                 xpb::Float64,ypb::Float64,zpb::Float64,xpe::Float64,ype::Float64,zpe::Float64,p::Param)

Input `kk' the vertical level, `fres' the Fresnel coefficient, `θ' the photon propagation direction, 
(`xpb',`ypb',`zpb') the coordinates of the photon, (`xpe',`ype',`zpe') the coordinates of the irradiance field,
and `p' the parameters.

Output the irradiance strength `ed1d', the information of the index `edi', `edj', and `edk', 
and the pointer `count', all as 1D arrays. 
"""
function energy!(ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},edj::Array{<:Int64,1},
                 edk::Array{<:Int64,1},count::Array{<:Int64,1},kk::Int64,fres::Float64,θ::Float64,
                 xpb::Float64,ypb::Float64,zpb::Float64,xpe::Float64,ype::Float64,zpe::Float64,p::Param)
    # check if downwelling
    if θ < pi/2
        # obtain crossing point on certain depth
        fract=(p.z[kk]-zpb)/(zpe-zpb)
        xpk=xpb+(xpe-xpb)*fract
        ypk=ypb+(ype-ypb)*fract    
        if xpk<p.xmin || xpk>=p.xmax || ypk<p.ymin || ypk>=p.ymax
            xpk,ypk=regulatexy(xpk,ypk,p)
        end
        
        # distribute ed to grid points
        area,i,j=interpolation2d(xpk,ypk,p)
        # distribute energy to grid point on level kk
        if count[1]==length(ed1d)
            error("Out of memory error! Increase the size of ed1d, edi, edj, and edk!")
        end
        for k=1:4
            "count is a dummy variable just to keep track of the size of the ed1d, edi, edj, edk"
            count[1]+=1
            "ed1d is the 1 by n matrix corresponding to the actual irradiance calculation"
            ed1d[count[1]]=(area[k]/p.dx/p.dy)*fres*cos(θ)
            "edi is the 1 by n matrix corresponding to the x coordination of the energy ex: ed(edi(count),edj(count),edk(count)) = ed1d(count)"
            edi[count[1]]=i[k]
            "edj is the 1 by n matrix corresponding to the y coordination of the energy ex: ed(edi(count),edj(count),edk(count)) = ed1d(count)"
            edj[count[1]]=j[k]
            "edk is the 1 by n matrix corresponding to the z coordination of the energy ex: ed(edi(count),edj(count),edk(count)) = ed1d(count)"
            edk[count[1]]=kk
        end
    end
    return nothing
end


function solar_energy!(esol::Array{<:Float64,2},kk::Int64,θ::Float64,
                       xpb::Float64,ypb::Float64,zpb::Float64,xpe::Float64,ype::Float64,zpe::Float64,
                       isign::Int64,p::Param)
    # initialize parameters
    dθ=pi/(p.num-1)
    # count the solar angle distribution of the radius
    if isign*(p.z[kk]-zpe) > 0
        dz1=abs(p.z[kk]-zpb)
        dz2=abs(p.z[kk]-zpe)
        w1=(floor(Int,dz1*abs(tan(θ))/p.dx)+1)/(p.x[p.nx]-p.x[1])
        w2=(floor(Int,dz2*abs(tan(θ))/p.dx)+1)/(p.x[p.nx]-p.x[1])
        if isign == 1            
            is=floor(Int,θ/dθ)+1
            res=mod(θ,dθ)/dθ
            if res > 0.5
                is=is+1       
            end
            esol[is,kk]=esol[is,kk]+w1
            esol[is,kk+1]=esol[is,kk+1]+w2
        end
        if isign == -1            
            is=floor(Int,θ/dθ)+1
            res=mod(θ,dθ)/dθ
            if res > 0.5
                is=is+1 
            end
            esol[is,kk+1]=esol[is,kk+1]+w1
            esol[is,kk]=esol[is,kk]+w2
        end
    else        
        w1=(floor(Int,abs(xpb-xpe)/p.dx)+1)/(p.x[p.nx]-p.x[1])
        if isign == 1            
            is=floor(Int,θ/dθ)+1
            res=mod(θ,dθ)/dθ
            if res > 0.5 
                is=is+1       
            end
            esol[is,kk]=esol[is,kk]+w1
        end          
        if isign == -1            
            is=floor(Int,θ/dθ)+1
            res=mod(θ,dθ)/dθ
            if res > 0.5
                is=is+1       
            end
            esol[is,kk+1]=esol[is,kk+1]+w1
        end      
    end
    return nothing
end

# function interaction(a::Float64,b::Float64,randrng)                 
#     # mu is cosine of scattering angle
#     ro1=rand(randrng)      
#     if ro1 > b/(a+b)         
#         θ=0.0
#         ϕ=0.0
#         isca=0
#         idie=1         
#     else         
#         ϕ=rand(randrng)*2*pi
#         ro2=rand(randrng)
#         mu=newton(ro2,mu)
#         θ=acos(mu)
#         isca=1
#         idie=0
#     end      
#     return θ,ϕ,isca,idie
# end

"""
    interaction(a, b, ph, θps, randrng)

Determine whether or not photon being absorb or scattering.
 
If the photon is scattered, determined the angle in the plane of the scattering event relative to the direction of photons before scattering
"""
function interaction(a::Float64,b::Float64,ph::Array{<:Float64,1},θps::Array{<:Float64,1},randrng)
    # this subroutine use petzold's scattering phase function
    ro1=rand(randrng)
    nums=size(ph,1)
    if ro1 > b/(a+b) 
        "single-scattering albedo, b/(a+b), determine the likeliness of being absorb or scatter. If ro1 is higher than single-scattering albedo, the photon is being absorbed"        
        θ=0.0
        "cosine of the scattering angle"
        ϕ=0.0
        isca=0
        idie=1         
    else         
        isca=1
        idie=0
        "azimuthal angle in the plane of the scattering event relative to the direction of photons before scattering"
        ϕ=rand(randrng)*2*pi
        ro2=rand(randrng)
        if ro2 < ph[1]
            #if the random number R is lower than the 0.517 (ph[1]), or the angle between 0 and 2.5 deg, (a little scattering)
            "treat the function of cumulative probabilty of scattering and corresponding angle (from 0 to 2.5 deg) as a linear functions"
            θ=(0+θps[1]*ro2/ph[1])
        else
            for i=1:nums
                #find the interval of cumulative probabilty of scattering (ph or ϕps) that would contain the random number R
                if ro2 > ph[i] && ro2 <= ph[i+1]
                    #the random number is in the interval (ph[i],ph[i+1]], ph[i] < R <= ph[i+1]
                    "interpolation function (?)"
                    θ=(θps[i]*(ph[i+1]-ro2)+θps[i+1]*(ro2-ph[i]))/(ph[i+1]-ph[i])
                    break
                end
            end
        end         
    end
    return θ,ϕ,isca,idie
end
      
# function newton(ro::Float64)      
#     #     use newton iteration to solve cosine of scattering angle
#     #     use pure water scattering phase function

#     f=0.835
#     itmax=10
#     a=-f/(6+2*f)
#     b=-3/(6+2*f)
#     c=0.5-ro    
#     mu=0
#     mu0=0
#     for it=1:itmax
#         fun=a*mu0^3+b*mu0+c
#         df=3*a*mu0^2+b
#         mu=mu0-fun/df
#         if abs(mu-mu0) < 1e-4
#             break
#         end
#         mu0=mu
#     end    
#     return mu
# end
        
function phaseFournier(μ::Float64,n::Float64)
    ν=(3-μ)/2
    ψ=[i.*0.25 for i=1:719].*π/180
    δ=(4/3/(n-1)^2).*(sin.(ψ/2)).^2
    δ180=(4/3/(n-1)^2)
    pf = (1 ./(4*π.*(1 .-δ).^2 .*δ.^ν).*(ν.*(1 .-δ)-(1 .-δ.^ν)+(δ.*(1 .-δ.^ν)-ν.*(1 .-δ))./(sin.(ψ/2)).^2) 
          .+ (1-δ180^ν).*(3 .*(cos.(ψ)).^2 .-1)./(16*π*(δ180-1)*δ180^ν))
    return pf2pfcum(pf,ψ),ψ[1:end-1]
end

function pf2pfcum(pf::Array{<:Float64,1},scag::Array{<:Float64,1})
    pfcum=cumsum(pf[1:end-1].*sin.(scag[1:end-1]).*diff(scag))./sum(pf[1:end-1].*sin.(scag[1:end-1]).*diff(scag))
    return pfcum
end

"""
    phasePetzold()

return 2 arrays: `ϕps` and `θps`. 

When `ϕps` is the cumulation distribution of scattering angle and
`θps` is the angle between new trajectory and the direction of the photon before scattering corresponding to each `ϕps`.
"""
function phasePetzold()
    #ph is the cumulative distribution function for the Petzold measurement
    # Kirk, 1981, Monte Carlo procedure for simulating the penetration of light into natural waters
    ϕps=reshape([0.517        0.666        0.751        0.811        0.851        0.879
        0.901        0.918        0.931        0.942        0.95        0.957
        0.963        0.968        0.972        0.975        0.978        0.981
        0.9832        0.9853        0.987        0.9888        0.9902        0.9917
        0.9928        0.9941        0.995        0.996        0.9968        0.9976
        0.9982        0.9988        0.9992        0.9997        0.9998        1.0]',(36,1))[:]          
    θps=reshape([2.5            7.5            12.5            17.5            22.5            27.5
            32.5            37.5            42.5            47.5            52.5            57.5
            62.5            67.5            72.5            77.5            82.5            87.5
            92.5            97.5            102.5            107.5            112.5            117.5
            122.5            127.5            132.5            137.5            142.5            147.5
            152.5            157.5            162.5            167.5            172.5            177.5]',(36,1))[:]   
    "Changing θps to radian unit" 
    θps = θps .* π/180
    return ϕps, θps
end

"""
    scatter(θ,ϕ,θs,ϕs) 

changing the scattered ray frame of referece, from the local system to the sphere coordination.
"""
function scatter(θ::Float64,ϕ::Float64,θs::Float64,ϕs::Float64)      
    "μx is the cartesian coordination in x direction of the unit vector in initial direction: (μx,μy,μz)"
    μx=sin(θ)*cos(ϕ)
    "μy is the cartesian coordination in y direction of the unit vector in initial direction: (μx,μy,μz)"
    μy=sin(θ)*sin(ϕ)
    "μz is the cartesian coordination in z direction of the unit vector in initial direction: (μx,μy,μz), when upward direction is positive"
    μz=-cos(θ)
    if abs(μz) > 0.99999  
        #Find cartesian coordination of the unit vector in new direction: (μxs,μys,μzs), if incoming light ray is in the downward direction
        "μxs is the cartesian coordination in x direction "
        μxs=sin(θs)*cos(ϕs)*μz/abs(μz)
        "μys is the cartesian coordination in y direction "
        μys=sin(θs)*sin(ϕs)*μz/abs(μz)
        "μzs is the cartesian coordination in z direction; when the upward direction is positive; last term would be postive if light coming upward"
        μzs=cos(θs)*μz/abs(μz)         
    else         
        μs=cos(θs)
        "sqms is the sin(θs)"
        sqms=(1-μs^2)^0.5
        "sqmz in the sin(θ)"
        sqmz=(1-μz^2)^0.5
        rsq=1/sqmz
        "changing the coordination from (θs,ϕs,r) in the local coordination to the (x,y,z) in cartesian coordination"
        μxs=μx*μz*rsq*sqms*cos(ϕs)-μy*rsq*sqms*sin(ϕs)+μx*μs
        μys=μy*μz*rsq*sqms*cos(ϕs)+μx*rsq*sqms*sin(ϕs)+μy*μs
        μzs=-sqmz*sqms*cos(ϕs)+μz*μs
    end
    "solar angle: angle between scattered ray and positive z axis(upward direction)"   
    θ=pi-acos(μzs)    
    μzs=(μxs^2+μys^2)^0.5
    "normalize the term μxs; -1 <= μxs <= 1; before normalize term μxs depend on z coordination of scattered ray"
    μxs=μxs/μzs
    "normalize the term μys; -1 <= μys <= 1; before normalize term μys depend on z coordination of scattered ray"
    μys=μys/μzs    
    if abs(μxs-1) < 1e-10
        #if μxs = 1, μys is approx 0, and the scattered ray is almost only in the direction of positive x 
        "azimuthal angle is zero, or the angle between light ray and positive x axis"
        ϕ=0.
    else
        if abs(μxs+1) < 1e-10
            #if μxs = -1, μys is approx 0, and the scattered ray is almost only in the direction of negative x 
            "azimuthal angle is 180 deg, or the angle between light ray and positive x axis"
            ϕ=pi
        else 
            if μys > 0 
                # if the scattered ray is in the first and second quadrant 
                "azimuthal angle is arccos of the x coordination"
                ϕ=acos(μxs)
            else 
                if μys < 0
                    # if the scattered ray is in the third and fourth quadrant 
                    "azimuthal angle is 360deg minus arccos of the x coordination"
                    ϕ=2*pi-acos(μxs)
                end
            end
        end
    end
    return Float64(θ),Float64(ϕ),Float64(θs),Float64(ϕs)
end

"""
    surface(idie::Int64,isca::Int64,xpe::Float64,ype::Float64,zpe::Float64,
                 area::Vector{Float64},interi::Vector{Int64},interj::Vector{Int64},
                 η::Array{<:Float64,2},ip::Int64,p::Param)
    
Checking whether or not the photon go through the surface back to air side
"""
function surface(idie::Int64,isca::Int64,xpe::Float64,ype::Float64,zpe::Float64,
                 area::Vector{Float64},interi::Vector{Int64},interj::Vector{Int64},
                 η::Array{<:Float64,2},ip::Int64,p::Param)

    # calculate surface elevation at (xpe,ype)
    interpolation2d!(area,interi,interj,xpe,ype,p)
    η0=0
    for k=1:4
        η0+=area[k]*η[interi[k],interj[k]]/p.dx/p.dy
    end
#    η0=(a1*η[i1,j1]+a2*η[i2,j2]+a3*η[i3,j3]+a4*η[i4,j4])/p.dx/p.dy
    # check if photon go back to air side
    if zpe >= η0
        idie=1
        isca=0
#        println("photon $ip go through surface back to air side!")
    else
        idie=0
    end    
    return idie,isca
end

"""
    interpolation2d(xpe, ype, p)

Calculate 4 value of areas between (xpe,ype), the coordination where photons end up, and different grid point around (xpe,ype)
"""
function interpolation2d(xpe::Float64,ype::Float64,p::Param)
    x=p.x
    y=p.y
    i=[0; 0; 0; 1]
    j=[0; 0; 0; 1]
    ""
    while x[i[4]] <= xpe && i[4] < size(x,1)
        "i[4] is the position in the grid in which the x-coordinate of our grid is larger than the x-coordination where photons end up with"
        i[4] += 1	
    end
    while y[j[4]] <= ype && j[4] < size(y,1)
        "j[4] is the position in the grid in which the y-coordinate of our grid is larger than the y-coordination where photons end up with"
        j[4] +=1
    end     
    i[1]=i[4]-1; j[1]=j[4]-1
    i[2]=i[1]+1; j[2]=j[1]     
    i[3]=i[1];   j[3]=j[1]+1
    "position of i,j that satisfied the boundary condition"
    i,j=regulateij(i,j,p.nx,p.ny)
    area=zeros(4)
    area[1]=(x[i[4]]-xpe)*(y[j[4]]-ype)
    area[2]=(xpe-x[i[3]])*(y[j[3]]-ype)
    area[3]=(x[i[2]]-xpe)*(ype-y[j[2]])
    area[4]=(xpe-x[i[1]])*(ype-y[j[1]])
    for k=1:4
        "in the case that area calculaion is less than 0, print the error"
        if area[k] < 0
            coorerror(area[k],i[k],j[k],xpe,ype,x,y,k)
        end
    end
    return area,i,j
end

function interpolation2d!(area::Vector{Float64},interi::Vector{Int64},
                          interj::Vector{Int64},xpe::Float64,ype::Float64,p::Param)
    x=p.x
    y=p.y
    interi[:]=[0; 0; 0; 1]
    interj[:]=[0; 0; 0; 1]
    while x[interi[4]] <= xpe && interi[4] < size(x,1)
        interi[4] += 1	
    end
    while y[interj[4]] <= ype && interj[4] < size(y,1)
        interj[4] +=1
    end     
    interi[1]=interi[4]-1; interj[1]=interj[4]-1
    interi[2]=interi[1]+1; interj[2]=interj[1]     
    interi[3]=interi[1];   interj[3]=interj[1]+1
    regulateij!(interi,interj,p.nx,p.ny)
    area[1]=(x[interi[4]]-xpe)*(y[interj[4]]-ype)
    area[2]=(xpe-x[interi[3]])*(y[interj[3]]-ype)
    area[3]=(x[interi[2]]-xpe)*(ype-y[interj[2]])
    area[4]=(xpe-x[interi[1]])*(ype-y[interj[1]])
    for k=1:4
        if area[k] < 0
            coorerror(area[k],interi[k],interj[k],xpe,ype,x,y,k)
        end
    end
    return nothing
end

"""
    coorerror(a ,i ,j ,xpe ,ype ,x ,y ,k)  

Produce an error, when our area surface calculation (from function interpolation2d) is lower than 0
"""
function coorerror(a::Float64,i::Int64,j::Int64,xpe::Float64,ype::Float64,x::Array{<:Float64,1},y::Array{<:Float64,1},k::Int64)
    println("In surface")
    println("a=$a")
    println("i=$i j=$j")
    println("x=$xpe y=$ype")
    println("x[i]=$(x[i]) y[j]=$(y[j])")
    println("k=$k")
    error("Coordinate calculation error!")
    return nothing
end

function updateed!(ed::Array{<:Float64,3},ed1d::Array{<:Float64,1},edi::Array{<:Int64,1},
                   edj::Array{<:Int64,1},edk::Array{<:Int64,1},n::Int64)
    for i=1:n
        ed[edi[i],edj[i],edk[i]] += ed1d[i]
    end
    return nothing
end

"""
    exported(ed::Array{<:Real,3},η::Array{<:Real,2},p::Param,
                  fname::String,mode="2D"::String,nk=0::Int64)

exported the irradiance data into `fname`.h5 file. 

There are 3 modes of exporting the data: `2D`, `3D`, and `full`. If not specified the mode will automatically set to `2D`. 
"""
function exported(ed::Array{<:Real,3},η::Array{<:Real,2},p::Param,
                  fname::String,mode="2D"::String,nk=0::Int64)
    if mode=="3D"
        # mode `3D` store only ed (3D irradiance grid) in the `fname` file
        h5open(fname*".h5","w") do fid
            write(fid,"ed",Float32.(ed))
        end
    end

    if mode=="full"
        # mode `full` store 4 dataset (the irradiance of the grid and 3 coordination array x,y,z) in the `fname` file
        x=similar(ed)
        y=similar(ed)
        z=similar(ed)
        for j=1:p.ny,k=1:p.nz
            # 3D array x stores the x coordination of the grid  
            x[:,j,k]=p.x        
        end
        for i=1:p.nx,k=1:p.nz
            # 3D array y stores the y coordination of the grid 
            y[i,:,k]=p.y
        end
        for i=1:p.nx,j=1:p.ny
            # 3D array z stores the z coordination of the grid 
            z[i,j,:]=p.z
        end
        h5open(fname*".h5","w") do fid
            write(fid,"ed",Float32.(ed))
            write(fid,"x",x)
            write(fid,"y",y)
            write(fid,"z",z)        
        end
    end

    "x2d stores the coordination of x at every depth"
    x2d=zeros(p.nx,p.nz)
    "y2d stores the coordination of y at every depth"
    y2d=zeros(p.ny,p.nz)
    "zx2d stores the z coordination (depth) at every x coordinattion"
    zx2d=zeros(p.nx,p.nz)
    "zy2d stores the z coordination (depth) at every y coordinattion"
    zy2d=zeros(p.ny,p.nz)
    for k=1:p.nz
        x2d[:,k]=p.x        
    end
    for k=1:p.nz
        y2d[:,k]=p.y
    end
    for i=1:p.nx
        zx2d[i,:]=p.z
    end
    for j=1:p.ny
        zy2d[j,:]=p.z
    end

    "xη is x coordination on the wave surface grid spanning from 0 to the total length of the whole grid "
    xη=[(i-1)*p.xl/p.nxη for i=1:p.nxs]
    "yη is y coordination on the wave surface grid spanning from 0 to the total length of the whole grid "
    yη=[(j-1)*p.yl/p.nyη for j=1:p.nys]
    h5open(fname*"xz.h5","w") do fid
        # `fname` + xz.h5 file stores 5 datasets: the irradiance field at 1st column, x2d, zx2d, xη and surface elevation at 1st column
        write(fid,"ed",ed[:,1,:])
        write(fid,"x",x2d)
        write(fid,"z",zx2d)
        write(fid,"xeta",xη)
        write(fid,"eta1D",η[:,1])
    end

    h5open(fname*"yz.h5","w") do fid
         # `fname` + yz.h5 file stores 5 datasets: the irradiance field at 1st row, y2d, zy2d, yη and surface elevation at 1st row
        write(fid,"ed",ed[1,:,:])
        write(fid,"y",y2d)
        write(fid,"z",zy2d)
        write(fid,"yeta",yη)
        write(fid,"eta1D",η[1,:])
    end

    if nk>0
        h5open(fname*"xy.h5","w") do fid
            # `fname` + xy.h5 file stores only the irradiance field at speific depth `nk`
            write(fid,"ed",ed[:,:,nk])
            # write(fid,"y",y2d)
            # write(fid,"z",zy2d)
            # write(fid,"yeta",yη)
            # write(fid,"eta1D",η[1,:])
        end
    end
    
    "the highest possible value of the number of grid that will be above the water surface"
    kcut=ceil(Int,(p.ztop-minimum(η))/p.dz)
    
    μ,σ,cv=edstats(ed,p)
    μ[1:kcut].=-1
    σ[1:kcut].=-1
    cv[1:kcut].=-1
    h5open(fname*"stats.h5","w") do fid
        #  `fname` + stats.h5 file stores 4 dataset: variance, mean, cv, and depth
        write(fid,"var",σ)
        write(fid,"mean",μ)
        write(fid,"cv",cv)
        write(fid,"z",p.z)
    end
    return nothing
end

"""
    edstats(ed::Array{<:Real,3},p::Param)

Calculate and return the value of mean `μ`, standard deviation `σ`, and `cv`, at each depth 
"""
function edstats(ed::Array{<:Real,3},p::Param)
    μ=zeros(p.nz)
    σ=zeros(p.nz)
    cv=zeros(p.nz)
    for k=1:p.nz
        σ2=mean(ed[:,:,k].^2)
        μ[k]=mean(ed[:,:,k])
        if μ[k]!=0
            cv[k]=sqrt(σ2/μ[k]^2-1)
        else
            cv[k]=-1
        end
        σ[k]=sqrt(σ2)
    end
    return μ,σ,cv
end

function ed2edr(ed::Array{<:Float64,3},nk::Int64,(iext,jext)::Tuple{<:Int64,<:Int64},p::Param)
    ed2d=ed[:,:,nk]./p.nphoton
    im=findmax(ed2d)[2][1]
    jm=findmax(ed2d)[2][2]
    edrm=maximum(ed2d)

    xind=(-iext:iext).*p.dx
    yind=(-jext:jext).*p.dy

    r=Float32[]
    edr=Float32[]
    regr=Float32.([i*p.dx for i=1:iext])
    regedr=zeros(Float32,iext)
    count=zeros(Int,iext)
    for i=1:iext*2+1, j=1:jext*2+1
        push!(r,sqrt(xind[i]^2+yind[j]^2))
        push!(edr,ed2d[im-iext-1+i,jm-jext-1+j])
        for k=2:iext-1
            if r[end]>=regr[k] && r[end]<regr[k+1]
                regedr[k]+=ed2d[im-iext-1+i,jm-jext-1+j]
                count[k]+=1
            end
        end
    end
    for k=1:iext
        if count[k]>0
            regedr[k]/=count[k]
        end
    end

    reged2d=copy(ed2d[im-iext:im+iext,jm-jext:jm+jext])
    for i=1:iext*2+1, j=1:jext*2+1
        tmpr=sqrt(xind[i]^2+yind[j]^2)
        if reged2d[i,j]<=1e-3*edrm
            for k=2:iext-1
                if tmpr>=regr[k] && tmpr<regr[k+1]
                    reged2d[i,j]=regedr[k]
                end
            end
        end
    end
    return im,jm,edrm,r,edr,ed2d[im-iext:im+iext,jm-jext:jm+jext],regr,regedr,reged2d
end