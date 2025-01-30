# Quick Start 

In this example, we introduce the `OceanLight` calculation of downwelling irradiance field underneath the flat surfcae. The code example here can be directly pasted onto Julia terminal, run through `.jl` script file, or IJulia notebook. 

First, we import `OceanLight` packages and another `random` dependent package, that will be used as a seed to generate random number in the Monte Carlo simulation.

```@setup  Center
cd(mktempdir()) 
using Pkg 
Pkg.add("OceanLight") 
Pkg.add("Plots") 
```
```@example Center 
using OceanLight 
using Random
```

##  Initial Condition 

OceanLight accesses all input variables in `.yml` format and stores their values in `Param` structure.

All input variables required by OceanLight can be separated into 3 categories:
1. **Irradiance:** Setting up the dimension of downwelling irradiance field calculation grid. 
2. **Photon:** Number of Photons and optical properties of the water.
3. **Wave:** Structure of surface wave field. 

```@example Center
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
```
To create a input variable file suitable for this package, user can either create a new file in `.yml` format, copy, and paste the code block above, or using a build-in function `OceanLight.writeparams` to automate the process. 

The function `OceanLight.writeparams` converts the dictionary of input variables into the `light.yml` file. 

```@example Center
data=Dict("irradiance"=>Dict("nxe"=>nxe,"nye"=>nye,"nz"=>nz,"dz"=>dz,"ztop"=>ztop,"num"=>num),
            "wave"=>Dict("pex"=>pex,"pey"=>pey,"nxeta"=>nxeta,"nyeta"=>nyeta),
            "photon"=>Dict("nxp"=>nxp,"nyp"=>nyp,"nphoton"=>nphoton,"a"=>a,"b"=>b,"kr"=>kr,"kbc"=>kbc))

OceanLight.writeparams(data)
```

Once the input variable file is generated, this package calculates other related variables and stores all values in the `Param` structure, through `OceanLight.readparams`. The default input variable file name that will be read is `light.yml`, but can be specified by users. In this case, we store `Param` structure in variable `p`. 

```@example Center
p = OceanLight.readparams()
```

## Initialize Parameters

Before the simulation, user needs to declare and initialize all parameters and their dimensions. 

### Air-Water Interaction

During the light refraction between two mediums calculation, `OceanLight` requires the surface elevation attribution: $\eta$; surface elevation, $\eta_{x}$; partial derivative of $\eta$ in x direction, and $\eta_{y}$; partial derivative of $\eta$ in y direction. All surface elevation distribution, as specified above, need to have the same dimension with the incoming photons' grid. Hence, 

```@example Center
η = zeros(p.nxs,p.nys)
ηx = zeros(p.nxs,p.nys)
ηy = zeros(p.nxs,p.nys)
```
After photons' interaction with the surface, OceanLight requires the information of specific coordinate of photon in cartesian grid ${xpb,ypb,zpb}$, the direction in which photon will travel in polar coordinate ${θ,ϕ}$, and the fraction of light that transmit through the water ${fres}$: all in the dimension of incoming photon grid size. 

```@example Center
xpb = zeros(p.nxp,p.nyp);
ypb = zeros(p.nxp,p.nyp);
zpb = zeros(p.nxp,p.nyp);
θ = zeros(p.nxp,p.nyp);
ϕ = zeros(p.nxp,p.nyp);
fres = zeros(p.nxp,p.nyp);
```

### Monte Carlo simulation

`OceanLight` tracks the path of each photon travelling inside the water medium and store the irradiance value in the grid `ed`. 

Users need to specify these variables and corresponding dimension. 

```@example Center
ed = zeros(p.nx, p.ny, p.nz)
esol = zeros(p.num, p.nz)
randrng = MersenneTwister(1234)
area = zeros(4)
interi = zeros(Int64,4)
interj = zeros(Int64,4)
ix = div(p.nxη,2)+1
iy = div(p.nyη,2)+1
ϕps,θps = OceanLight.phasePetzold()
```

## Air-Water Interaction

During the air-water interaction process, OceanLight simulates the photons transfer directly downward from the air side, interacts with the water surface, and transfer down into water medium. 

User can generate random surface elevation attribution ${\eta,\eta_{x},\eta_{y}}$ with `OceanLight.setwave!`, or provided specific data ${\eta_{0},\eta_{x0},\eta_{y0}}$  . OceanLight can map the user's provided data of ${\eta_{0},\eta_{x0},\eta_{y0}}$, which might have different dimension onto the suitable dimension of input value ${η,ηx,ηy}$ with `OceanLight.convertwave!`.

In this example, we will consider the case of flat surface elevation. Hence, ${η,ηx,ηy}$ is equal to the matrix of zeros. 

Once all the input variables are in place, `OceanLight.interface!` calculate the refraction of the light between two medium given surface elevation attribution and return the position, reflectance angle, and transmission ratio. 

```@example Center
OceanLight.interface!(xpb,ypb,zpb,θ,ϕ,fres,η,ηx,ηy,p)
```

## Monte Carlo simulation

`OceanLight` simulates the photon traveling inside the water medium, given its initial position ${xpb,ypb,zpb}$ and the direction it started with${θ,ϕ}$. Once photons are inside the water, `OceanLight` will track its path, governed by its probability distribution and the attenuated coefficient input, and store the irradiance value in the grid `ed`. 

The `transfer!` function simulate a single photon path and store its landed position on the grid `ed`. Hence, to simulate multiple photons, users need to loop the `transfer!` function and giving the input of an individual photon's number `ip`. Thus, `OceanLight` could facilitate parallel computation. 

```@example Center
for ip = 1:p.nphoton
    OceanLight.transfer!(ed,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,xpb[ix,iy],
        ypb[ix,iy],zpb[ix,iy],area,interi,interj,randrng,η,ϕps,θps,p,1)
end
```

Lastly, once the field `ed` is obtained. The `applybc!` apply and ensure the boundary condition of the `ed`. 

```@example Center
OceanLight.applybc!(ed,p)
```

Once the simulation is over, the downwelling irradiance field $I(x,y,z)$ is stored in variable `ed`. To export the data, the build-in function `OceanLight.exported` export data and its statistical information in to `.h5` format, depended on the path given by user. 

The result downwelling irradiance field $I(x,y,z)$ is in 3 dimension tensor, where the first two dimensions represent horizontal field, and the last represents the vertical field or depth. The first few depth layer represents the air phase, hence $I(x,y,z)$ is zeros, depending on how many layer specified as air phases `ztop`. To visualize using `Plots` julia package, users can try the example code below. 

```@example Center
using Plots 
using Plots.Measures

l = @layout [grid(2,1) a{0.5w} ; b{0.5w}]

p1 = heatmap(p.x,p.y,log.(ed[:,:,40]),clim=(-20,20),framestyle = :box,grid = false, c =cgrad(:viridis)
    ,legend = :none,xlabel="\$x(m)\$",ylabel="\$y(m)\$")
p2 = heatmap(p.x,p.y,log.(ed[:,:,160]),clim=(-20,20),framestyle = :box,grid = false, c =cgrad(:viridis)
    ,legend = :none,xlabel="\$x(m)\$",ylabel="\$y(m)\$")
p3 = heatmap(p.x,p.z,reverse(transpose(log.(ed[:,256,:]))),clim=(-20,20),framestyle = :box,grid = false, c =cgrad(:viridis)
    ,xlabel="\$x(m)\$",ylabel="\$z(m)\$"
    ;cbar_title="\$\\ln(I(x,y,z))\$")
plot(p1, p2,p3, layout = l,
title = ["($i)" for j in 1:1, i in ["a","b","c"]], titleloc = :left, titlefont = font(8)
,left_margin = [10mm 0mm],right_margin = [10mm 0mm])
plot!(size=(900,1200))
```

