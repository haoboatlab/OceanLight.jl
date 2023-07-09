using SerialFFT
using Random
using Statistics
using HDF5
include("LightMC.jl")
include("spec.jl")
using MPI
MPI.Init()
comm = MPI.COMM_WORLD
myid = MPI.Comm_rank(comm)
ncpu = MPI.Comm_size(comm)

p=readparams()

println("This program generates the raw data of a single photon from multiple surface wave data")

nind=p.nphoton
if mod(nind,ncpu) ==0
    dind=div(nind,ncpu)
else
    dind=div(nind,ncpu) + 1
end
inds=myid*dind+1
inde=(myid+1)*dind
if inde > nind
    inde=nind
end
#allind=CartesianIndices(zeros(p.nphoton))
allind=1:p.nphoton
η=zeros(p.nxs,p.nys)
ηx=zeros(p.nxs,p.nys)
ηy=zeros(p.nxs,p.nys)
randrng = MersenneTwister(1234+myid)
ϕps,θps=phasePetzold()

ed=zeros(p.nx,p.ny,p.nz)
esol=zeros(p.num,p.nz)
area=zeros(4)
interi=zeros(Int64,4)
interj=zeros(Int64,4)
xpb=zeros(p.nxp,p.nyp)
ypb=zeros(p.nxp,p.nyp)
zpb=zeros(p.nxp,p.nyp)
θ=zeros(p.nxp,p.nyp)
ϕ=zeros(p.nxp,p.nyp)
fres=zeros(p.nxp,p.nyp)
ix=div(p.nxη,2)+1
iy=div(p.nyη,2)+1

icase1=parse(Int,ARGS[1])
icase2=parse(Int,ARGS[2])

for icase=icase1:icase2
    η.=0
    ηx.=0
    ηy.=0
    ed.=0
    esol.=0    
    if myid==0
        run(`mkdir -p rawdat/case$(icase)`)
        fid=h5open("allparams/surfwave$(icase).h5","r")
        η0=read(fid,"eta")
        ηx0=read(fid,"ex")
        ηy0=read(fid,"ey")
        close(fid)
        convertwave!(η,ηx,ηy,η0,ηx0,ηy0,p.kbc)
    end
    MPI.Allreduce!(η,+,comm)
    MPI.Allreduce!(ηx,+,comm)
    MPI.Allreduce!(ηy,+,comm)
    if p.z[1] <= maximum(η)
        error("ztop smaller than maximum η!")
    end
    println("myid $myid wave surface input complete!")
    
    interface!(xpb,ypb,zpb,θ,ϕ,fres,η,ηx,ηy,p)
    println("myid $myid refraction at interface complete!")
    @time begin
        for ind=inds:inde
            ip=allind[ind]
            transfer!(ed,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,
                      xpb[ix,iy],ypb[ix,iy],zpb[ix,iy],area,interi,interj,
                      randrng,η,ϕps,θps,p,1)
        end
    end

    MPI.Allreduce!(ed,+,comm)
    
    if myid ==0
        applybc!(ed,p)
        exported(ed,η,p,"rawdat/case$(icase)/ed","3D",440)
    end
    println("$myid completes case $icase !")
end

println("$myid Program ends!")

