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
comm = MPI.COMM_WORLD
# myid is rank of each cpu (tagname corresponding to each cpu being used)
myid = MPI.Comm_rank(comm)
# ncpu is number of cpu used to run the program
ncpu = MPI.Comm_size(comm)
```
```@example MPI_tutorial 
# the total number of photons in this case is the photons at each grid poin multiply by every single grid
nind=p.nphoton
# dind is the number photons that each cpu need to work on

if mod(nind,ncpu) ==0
    dind=div(nind,ncpu)
else
    dind=div(nind,ncpu) + 1
end
inds=myid*dind+1
# inde is the total number of the photons that cpu would simulate
inde=(myid+1)*dind
# if total number is larger than the actual photons we want to calculate, set it equal

if inde > nind
    inde=nind
end

allind=1:p.nphoton
"η is the height(z axis) corresponding to each grid point on the surface, 2d array of grid number in \
x and y direction"
η=zeros(p.nxs,p.nys)
"ηx is the x coordination corresponding to each grid point on the surface, 2d array of grid number in\
 x and y direction"
ηx=zeros(p.nxs,p.nys)
"ηy is the y coordination corresponding to each grid point on the surface, 2d array of grid number in\
 x and y direction"
ηy=zeros(p.nxs,p.nys)
randrng = MersenneTwister(1234+myid)
ϕps,θps=phasePetzold()

ed=zeros(p.nx,p.ny,p.nz)
esol=zeros(p.num,p.nz)
area=zeros(4)
interi=zeros(Int64,4)
interj=zeros(Int64,4)
ix=div(p.nxη,2)+1
iy=div(p.nyη,2)+1


η.=0
ηx.=0
ηy.=0
ed.=0
esol.=0
MPI.Allreduce!(η,+,comm)
MPI.Allreduce!(ηx,+,comm)
MPI.Allreduce!(ηy,+,comm)
if p.z[1] <= maximum(η)
    error("ztop smaller than maximum η!")
end
println("myid $myid wave surface input complete!")
xpb,ypb,zpb,θ,ϕ,fres=interface(η,ηx,ηy,p)
println("myid $myid refraction at interface complete!")

for ind=inds:inde
    ip=allind[ind]
    transfer!(ed,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,
                xpb[ix,iy],ypb[ix,iy],zpb[ix,iy],area,interi,interj,
                                    randrng,η,ϕps,θps,p,1)
end


MPI.Allreduce!(ed,+,comm)
```

```
if myid ==0
    applybc!(ed,p)
    exported(ed,η,p,"rawdat/case$(icase)/ed","3D",176)
end
```
```@example MPI_tutorial 
MPI.Barrier(comm)
MPI.Finalize()
```

# Run .jl file

```
mpirun -np 4 julia FILE_NAME.jl
```

# Visualization 







```@example MPI_tutorial 
using Plots
using Plots.Measures

# Choose slice indices
iy_c = 256   # y-index for vertical cross-section
iz_a = 40    # z-index for panel (a)
iz_b = 160   # z-index for panel (b)

# Define layout: 2 rows, 2 columns, but right column spans both rows
l = @layout [grid(2,1) c]

# Panel (a) : z = iz_a
p1 = heatmap(
    p.x .-10, p.y .-10, log.(ed[:,:,iz_a]),
    clim=(-20,-7), framestyle=:box, grid=false,
    c=cgrad(:viridis), legend=:none,
    xlabel="\$x(m)\$", ylabel="\$y(m)\$",
    title="(a) z = $(round(p.z[iz_a], digits=1)) m",
    xlim=[minimum(p.x).-10, maximum(p.x).-10],
    ylim=[minimum(p.y).-10, maximum(p.y).-10])
plot!(p1, [minimum(p.x).-10, maximum(p.x).-10], [p.y[iy_c]-10, p.y[iy_c]-10],
      color=:red, lw=2, ls=:dash, alpha=0.6)

# Panel (b) : z = iz_b
p2 = heatmap(
    p.x .-10, p.y .-10, log.(ed[:,:,iz_b]),
    clim=(-20,-7), framestyle=:box, grid=false,
    c=cgrad(:viridis), legend=:none,
    xlabel="\$x(m)\$", ylabel="\$y(m)\$",
    title="(b) z = $(round(p.z[iz_b], digits=1)) m",
    xlim=[minimum(p.x).-10, maximum(p.x).-10],
    ylim=[minimum(p.y).-10, maximum(p.y).-10])
plot!(p2, [minimum(p.x).-10, maximum(p.x).-10], [p.y[iy_c]-10, p.y[iy_c]-10],
      color=:red, lw=2, ls=:dash, alpha=0.6)

# Panel (c) : vertical cross-section at y = iy_c
p3 = heatmap(
    p.x .-10, reverse(p.z), reverse(transpose(log.(ed[:,iy_c,:]))),
    clim=(-20,-7), framestyle=:box, grid=false,
    c=cgrad(:viridis), ylim=(-(nz*dz-10),0),
    xlabel="\$x(m)\$", ylabel="\$z(m)\$",
    cbar_title="\$\\ln\\frac{I(x,y,z)}{I_{0}}\$",
    title="(c) y = $(round(p.y[iy_c]-10, digits=1)) m",
    xlim=[minimum(p.x).-10, maximum(p.x).-10])

# Add horizontal lines for z = iz_a, iz_b
plot!(p3, [minimum(p.x).-10, maximum(p.x).-10], [p.z[iz_a], p.z[iz_a]],
      color=:red, lw=1.5, ls=:dash, label="", alpha=0.6)
plot!(p3, [minimum(p.x).-10, maximum(p.x).-10], [p.z[iz_b], p.z[iz_b]],
      color=:red, lw=1.5, ls=:dash, label="", alpha=0.6)

# Combine
plot(p1, p2, p3, layout=l,
     size=(900,700),
     titleloc=:left, titlefont=font(8),
     left_margin=10mm, right_margin=10mm)

```