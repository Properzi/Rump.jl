using Rump
using Test

@testset "greet" begin
    @test greet_Rump() == "Hello Rump!"
    @test greet_Rump() != "Hello world!"
end

@testset "LAlgebra" begin
    M=[2 2; 1 2]
    N=[1 1; 1 2]
    P=[3 1 3; 3 3 3; 1 2 3]
    @test check_LAlgebra(M) == true
    @test check_LAlgebra(N) == false
    #@test LAlgebra(M, check=true) == LAlgebra([2 2; 1 2])
    @test LAlgebra(M) == LAlgebra(M)
    @test_throws ErrorException("[1 1; 1 2] does not define an L-algebra") LAlgebra(N, check=true)
    
    A=LAlgebra(M)
    B=LAlgebra(P)
    x=LAlgebraElem(A,1)
    y=LAlgebraElem(A,2)
    b=LAlgebraElem(B,2)
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