using LightMC
using Test 
using YAML
using HDF5
using Random
using UnicodePlots
using Statistics
using MPI 

MPI.Init()
comm = MPI.COMM_WORLD
myid = MPI.Comm_rank(comm)
ncpu = MPI.Comm_size(comm)

@testset "initialize MPI" begin
    @test myid == MPI.Comm_rank(comm)
    @test ncpu == MPI.Comm_size(comm)
end

p = LightMC.readparams("data/initial_condition/multipleCPU/light.yml")
ϕps,θps = LightMC.phasePetzold()
allind=CartesianIndices((1:p.nxp,1:p.nyp,1:p.nphoton))

η=zeros(p.nxs,p.nys)
ηx=zeros(p.nxs,p.nys)
ηy=zeros(p.nxs,p.nys)
seedid=rand(1:2000,1)[1]
randrng = MersenneTwister(seedid+myid)  

if myid==0
    fid=h5open("data/initial_condition/multipleCPU/surfwave.h5","r")
    η0=read(fid,"eta")
    ηx0=read(fid,"ex")
    ηy0=read(fid,"ey")
    close(fid)
    convertwave!(η,ηx,ηy,η0,ηx0,ηy0,p.kbc)  
end
MPI.Allreduce!(η,+,comm)
MPI.Allreduce!(ηx,+,comm)
MPI.Allreduce!(ηy,+,comm)

nind=p.nxp*p.nyp*p.nphoton
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

@testset "Initial condition for MPI" begin
    @test dind >= div(nind,ncpu)
    @test inds > 0 
end

if myid==0
    ed=zeros(p.nx,p.ny,p.nz)
end
esol=zeros(p.num,p.nz)

ncont=10000000
ed1d=zeros(Float64,ncont)
edi=zeros(Int,ncont)
edj=zeros(Int,ncont)
edk=zeros(Int,ncont)
count=[0]

nout=10000
iout=0
xpb,ypb,zpb,θ,ϕ,fres=interface(η,ηx,ηy,p)

@time begin
    for ind=inds:inde
        ix=allind[ind][1]
        iy=allind[ind][2]
        ip=allind[ind][3]
        transfer!(ed1d,edi,edj,edk,count,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,
                  xpb[ix,iy],ypb[ix,iy],zpb[ix,iy],randrng,η,ϕps,θps,p,1)

        if mod(ind-inds+1,nout)==0 || count[1]>=ncont/2
            if myid!=0
                MPI.Send(count,0,myid,comm)
                if count[1]>0
                    MPI.Send(ed1d[1:count[1]],0,10000+myid,comm)
                    MPI.Send(edi[1:count[1]],0,20000+myid,comm)
                    MPI.Send(edj[1:count[1]],0,30000+myid,comm)
                    MPI.Send(edk[1:count[1]],0,40000+myid,comm)
                    count[1]=0
                end
            else
                updateed!(ed,ed1d,edi,edj,edk,count[1])
                count[1]=0
                for i=1:ncpu-1
                    ndat=[0]
                    MPI.Recv!(ndat,i,i,comm)
                    if ndat[1]>0
                        MPI.Recv!(ed1d,i,10000+i,comm)
                        MPI.Recv!(edi,i,20000+i,comm)
                        MPI.Recv!(edj,i,30000+i,comm)
                        MPI.Recv!(edk,i,40000+i,comm)
                        updateed!(ed,ed1d,edi,edj,edk,ndat[1])
                    end
                end
            end
        end
    end

    if myid!=0
        MPI.Send(count,0,myid,comm)
        if count[1]>0
            MPI.Send(ed1d[1:count[1]],0,10000+myid,comm)
            MPI.Send(edi[1:count[1]],0,20000+myid,comm)
            MPI.Send(edj[1:count[1]],0,30000+myid,comm)
            MPI.Send(edk[1:count[1]],0,40000+myid,comm)
            count[1]=0
        end
    else
        updateed!(ed,ed1d,edi,edj,edk,count[1])
        count[1]=0
        for i=1:ncpu-1
            ndat=[0]
            MPI.Recv!(ndat,i,i,comm)
            if ndat[1]>0
                MPI.Recv!(ed1d,i,10000+i,comm)
                MPI.Recv!(edi,i,20000+i,comm)
                MPI.Recv!(edj,i,30000+i,comm)
                MPI.Recv!(edk,i,40000+i,comm)
                updateed!(ed,ed1d,edi,edj,edk,ndat[1])
            end
        end
    end
end

test_data = tempdir()
applybc!(ed,p)
exported(ed,η,p,test_data*"/ed","3D")

Test_ed = h5open(test_data*"/ed.h5","r")
test_ed = read(Test_ed,"ed")
close(Test_ed)

Test_edstats = h5open(test_data*"/edstats.h5","r")
test_cv = read(Test_edstats, "cv")
test_mean = read(Test_edstats, "mean")
test_var = read(Test_edstats, "var")
test_z = read(Test_edstats, "z")
close(Test_edstats)

Test_edxz = h5open(test_data*"/edxz.h5","r")
test_edxz = read(Test_edxz,"ed")
close(Test_edxz)

Test_edyz = h5open(test_data*"/edyz.h5","r")
test_edyz = read(Test_edyz,"ed")
close(Test_edyz)

Bench_mean = h5open("data/benchmark_data/multipleCPU/ed_mean.h5","r")
bench_mean = read(Bench_mean,"mean")
close(Bench_mean)

MaxPercentDif(array1,array2) = findmax(broadcast(abs,((array1-array2)/(array1))*100))
(Diff_mean,loc_mean) = MaxPercentDif(test_mean,bench_mean)

@testset "Result" begin
    @testset "export data" begin
        @test length(test_z) == parameter.nz
        @test mean(ed[:,:,50]) == test_mean[50]
        @test sqrt(mean(ed[:,:,30].^2)) == test_var[30]
        @test test_cv[1] == -1 
        @test floor(sqrt(test_var[60]^2/test_mean[60]^2-1)) == floor(test_cv[60])
        @test last(test_ed[1,:,:]) == last(test_edxz[:,:])
        @test last(test_ed[:,1,:]) == last(test_edyz[:,:])
    end 
    @testset "comparison with benchmark" begin
        @test Diff_mean <= 2
        @test abs((mean(ed[:,:,50])-bench_mean[50])/(bench_mean[50])*100) <= 2
    end 
end