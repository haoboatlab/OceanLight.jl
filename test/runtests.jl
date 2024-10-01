using LightMC
using Test 

@testset "LightMC.jl" begin
    @info "testing the functions being used in lightMC.jl package"
    run(`julia test_function.jl`)
end
