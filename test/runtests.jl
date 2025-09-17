using Rump
using Test

@testset "greet" begin
    @test greet_Rump() == "Hello Rump!"
    @test greet_Rump() != "Hello world!"
end

@testset "l_algebra" begin
    M=[2 2; 1 2]
    N=[1 1; 1 2]
    P=[3 1 3; 3 3 3; 1 2 3]
    @test check_l_algebra(M) == true
    @test check_l_algebra(N) == false
    #@test l_algebra(M, check=true) == l_algebra([2 2; 1 2])
    @test l_algebra(M) == l_algebra(M)
    @test_throws ErrorException("[1 1; 1 2] does not define an L-algebra") l_algebra(N, check=true)
    
    A=l_algebra(M)
    B=l_algebra(P)
    x=l_algebra_element(A,1)
    y=l_algebra_element(A,2)
    b=l_algebra_element(B,2)
    @test x*x == y
    @test size(A) == 2
    @test length(A) == 2
    @test x in A && b in B && !(b in A)
    
    e = logical_unit(A)
    @test e == y
    @test x*x == x* e && x*e == e && e * x == x
    #@test 
    #@test 
    #@test 
    #@test 
end