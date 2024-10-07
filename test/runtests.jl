using OceanLight
using Test 

@testset "OceanLight.jl" begin
    @info "testing the functions being used in OceanLight.jl package"
    run(`julia test_function.jl`)
end
