# Implementing MPI

Since `OceanLight.jl` implements the Monte Carlo ray tracing algorithm, the accuracy of the solution depends on its convergence. However, this process can be inefficient and time-consuming for certain problems that require a large number of photons `nphoton`. To address this issue, instead of running the entire simulation on a single CPU, we can divide the workload by allocating smaller chunks of photons across multiple CPUs, simulate them independently, and then combine the results into a single solution. That’s where MPI comes in.

MPI, or Message Passing Interface, is a communication library designed for sending and receiving messages across multiple processors, making it suitable for parallel programming. To install MPI on your system, check out the [Installing Open MPI guide](https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html)

```@setup  MPI_tutorial
cd(mktempdir()) 
using Pkg 
Pkg.add("OceanLight") 
Pkg.add("Plots") 
Pkg.add("MPI")
```

## Problem

In this example, the goal is to calculate the downwelling irradiance field when the surface is completely flat, and a total of 10,000,000 photons are focused at a single point at the center. The domain of interest is defined as $x,y \in \left[\mathrm{-10m},\mathrm{10m}\right]$, and $z \in \left[\mathrm{-190m},\mathrm{10m}\right]$ in depth, corresponding to a grid resolution of $512 \times 512 \times 200$ points. Periodic boundary conditions are applied at the domain boundaries. The attenuation properties of water are characterized by an absorption coefficient of $a = 0.0196$ and a scattering coefficient of $b = 0.0031$, which correspond to the optical properties of seawater at a wavelength of $490 ,\mathrm{nm}$ [^1].

To execute parallel programming using MPI, the Julia program must be run from a `.jl` script file. Inside the script, we first create the input data structures for the simulation. Here, `OceanLight.writeparams` generates a new input file in `.yml` format.

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

For a step-by-step description and the utility of each function, check out this [tutorial](https://haoboatlab.github.io/OceanLight.jl/dev/QuickStart/Center/). 

## Initialize MPI and variables

To enable parallel programming in Julia, we first import the `MPI` package and then initialize the environment with `MPI.Init()`.  

```@example MPI_tutorial 
using MPI
MPI.Init()
```

Here, `MPI.COMM_WORLD` is a predefined intracommunicator that represents the group of processes launched with MPI. We assign it to the variable `comm`. Each processor within the communicator has a unique rank (or ID number), which can be accessed using `MPI.Comm_rank(comm)`. Finally, the total number of processors in the communicator can be obtained with `MPI.Comm_size(comm)`. 

```@example MPI_tutorial 
comm = MPI.COMM_WORLD
processor_id = MPI.Comm_rank(comm)
number_cpu = MPI.Comm_size(comm)
```

Users need to specify these variables and their corresponding dimensions.

**NOTE:** The random number generator is created with `MersenneTwister`. Given the same input or seed, this pseudorandom number generator (PRNG) will always produce random numbers in a predictable and reproducible way. Hence, to ensure the independence of Monte Carlo paths across all processors, we need to offset the seed in each processor so that it does not generate the exact same random numbers.`ix` and `iy` are the array coordination corresponding to the center of the field.  

```@example MPI_tutorial 
p = OceanLight.readparams()

η = zeros(p.nxs,p.nys)
ηx = zeros(p.nxs,p.nys)
ηy = zeros(p.nxs,p.nys)
randrng = MersenneTwister(1234 + processor_id)
ϕps,θps = phasePetzold()
ed = zeros(p.nx,p.ny,p.nz)
esol = zeros(p.num,p.nz)
area = zeros(4)
interi = zeros(Int64,4)
interj = zeros(Int64,4)
ix = div(p.nxη,2) + 1
iy = div(p.nyη,2) + 1
```

## Monte Carlo path with MPI

The general idea behind implementing parallel programming in `OceanLight.jl` is that multiple CPUs can operate simultaneously to obtain results. First, we define a range of photon IDs, `photon_list`, which is a sequence of positive integers representing the total number of photons to be simulated. From this range, an equal number of photons are distributed to each processor and stored as `cpu_number_photon`. Once the refraction at the interface has been calculated, all photons proceed to follow independent Monte Carlo paths. We then iterate through each photon on every processor, specifying its starting point `starting` and ending point `ending` within the `for` loop. Finally, after all simulations are completed, the results from all processors are combined into a single output. 

```@example MPI_tutorial 
photon_list = 1:p.nphoton
total_number_photon = p.nphoton

if mod(total_number_photon,number_cpu) ==0
    cpu_number_photon = div(total_number_photon , number_cpu)
else
    cpu_number_photon = div(total_number_photon , number_cpu) + 1
end

starting = (processor_id * cpu_number_photon) + 1
ending = (processor_id + 1) * cpu_number_photon

if ending > total_number_photon
    ending = total_number_photon
end

MPI.Allreduce!(η,+,comm)
MPI.Allreduce!(ηx,+,comm)
MPI.Allreduce!(ηy,+,comm)
if p.z[1] <= maximum(η)
    error("ztop smaller than maximum η!")
end

xpb,ypb,zpb,θ,ϕ,fres = interface(η,ηx,ηy,p)

for ind = starting:ending
    ip = photon_list[ind]
    OceanLight.transfer!(ed,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,
                xpb[ix,iy],ypb[ix,iy],zpb[ix,iy],area,interi,interj,
                                    randrng,η,ϕps,θps,p,1)
end

println("Processor ID $processor_id complete!")
MPI.Allreduce!(ed,+,comm)
```

We apply periodic boundary conditions and export the results in `.h5` format for later accesability. 

```
if processor_id == 0
    OceanLight.applybc!(ed,p)
    OceanLight.exported(ed,η,p,"ed","3D",176)
end
```
Once `MPI` finishes its task, we end the communication process with `MPI.Finalize`. 

```@example MPI_tutorial 
MPI.Barrier(comm)
MPI.Finalize()
OceanLight.applybc!(ed,p) # hide
```

# Run .jl file

Once the Julia script file is saved, the code is ready to be executed. Within the same folder, use the command line to run `mpirun` and execute the MPI task. Here, the number following `-np` icorresponds to the number of CPUs you want to use, and do not forget to replace `FILE_NAME.jl` with the name of your Julia script file.

```
mpirun -np 4 julia FILE_NAME.jl
```

# Visualization 

To visualize the downwelling irradiance field, we first extract the data from the `ed.h5` file, along with the spatial coordinates in the x, y, and z directions, into the `struct p`. The code example below can be run directly in the Julia terminal, executed from a `.jl` script file, or used within an IJulia notebook.

```
using OceanLight 
using HDF5

p = OceanLight.readparams()
ed = h5open("ed.h5" ,"r")
ed = read(ed,"ed")
```

The visualization process is similar to the example shown [here](https://haoboatlab.github.io/OceanLight.jl/dev/QuickStart/Center/). 

```@example MPI_tutorial 
using Plots
using Plots.Measures

max_val, max_loc = findmax(ed)
ed = ed./max_val
nonzero_vals = ed[ed .!= 0]
min_val = minimum(nonzero_vals)

for i in 1:Int(nxe+1)
    for j in 1:Int(nye+1)
        for k in ztop:nz
            if ed[i,j,k] == 0
                ed[i,j,k] = min_val
            end
        end
    end
end

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

## Reference 

[^1]: Smith, R. C., & Baker, K. S. (1981). Optical properties of the clearest natural waters (200-800 nm). Applied optics, 20(2), 177–184. https://doi.org/10.1364/AO.20.000177 