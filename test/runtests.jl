using Rump
using Test

@testset "greet" begin
    @test Rump.greet_Rump() == "Hello Rump!"
    @test Rump.greet_Rump() != "Hello world!"
end

@testset "check" begin
    @test Rump.check_LAlgebra([2 2; 1 2]) == true
    @test Rump.check_LAlgebra([1 2; 1 2]) == false
end