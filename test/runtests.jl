using LightMC
using Test 
using YAML
using MPI
using UnicodePlots

@testset "LightMC.jl" begin
    @info "testing the functions being used in lightMC.jl package"
    run(`julia test_function.jl`)
    
    @info "testing muliple cpu via MPI package"
    run(`mpirun -n 4 julia test_mpi.jl`)
end
