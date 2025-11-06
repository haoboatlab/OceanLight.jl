using HydrOptics
using Test 
using YAML
using HDF5
using Random
using Statistics

parameter = HydrOptics.readparams("data/initial_condition/singleCPU/light.yml")
ϕps,θps = HydrOptics.phasePetzold()

η=zeros(parameter.nxs,parameter.nys)
ηx=zeros(parameter.nxs,parameter.nys)
ηy=zeros(parameter.nxs,parameter.nys)

fid=h5open("data/initial_condition/singleCPU/surfwave1.h5","r")
η0=read(fid,"eta")
ηx0=read(fid,"ex")
ηy0=read(fid,"ey")
close(fid)
    
HydrOptics.convertwave!(η,ηx,ηy,η0,ηx0,ηy0,parameter.kbc)

@testset "Initial Condition" begin
    @testset "readparams()" begin
        @test parameter.kbc >= 0
    end
    @testset "convertwave()" begin
        if parameter.kbc == 0
            @test ηy[end,end] == ηy[1,1]
        else
            @test ηy[1:size(η0,1),1:size(η0,2)] == ηy0
        end
    end
    @testset "phasePetzold()" begin
        @test size(ϕps) == (36,)
        @test size(θps) == (36,)
    end
end

xpb,ypb,zpb,θ,ϕ,fres = HydrOptics.interface(η,ηx,ηy,parameter)

@testset "Refraction between atmosphere and water" begin
    @testset "interface()" begin
        @test size(fres) == (parameter.nxp,parameter.nyp)
        @test size(θ) == (parameter.nxp,parameter.nyp)
        @test size(ϕ) == (parameter.nxp,parameter.nyp)
    end 
end

allind = 1:parameter.nphoton
ed = zeros(parameter.nx, parameter.ny, parameter.nz)
esol = zeros(parameter.num, parameter.nz)
randrng = MersenneTwister(1234)
area=zeros(4)
interi=zeros(Int64,4)
interj=zeros(Int64,4)
ix=div(parameter.nxη,2)+1
iy=div(parameter.nyη,2)+1

@time begin
    for ind=1:parameter.nphoton
        ip=allind[ind]
        HydrOptics.transfer!(ed,esol,θ[ix,iy],ϕ[ix,iy],fres[ix,iy],ip,xpb[ix,iy],
        ypb[ix,iy],zpb[ix,iy],area,interi,interj,randrng,η,ϕps,θps,parameter,1)
    end
end

@testset "Monte Carlo Simulation" begin
    @testset "transfer()" begin
        @test ed[ix,iy,floor(Int, parameter.ztop)+2] >= 0
        @test ed[ix,iy,floor(Int, parameter.ztop)+3] >= 0
        @test ed[ix,iy,floor(Int, parameter.ztop)+4] >= 0
        @test ed[ix,iy,floor(Int, parameter.ztop)+5] >= 0
    end 
end

test_data = tempdir()
HydrOptics.applybc!(ed,parameter)
HydrOptics.exported(ed,η,parameter,test_data*"/ed","3D")

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

Bench_mean = h5open("data/benchmark_data/singleCPU/ed_mean.h5","r")
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
