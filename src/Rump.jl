@doc """
Welcome to Rump version $(VERSION_NUMBER)

Rump is developed...bla bla

"""
module Rump


include("imports.jl")
const VERSION_NUMBER = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "Project.toml"))["version"])

export greet_Rump
export check_LAlgebra
include("functions.jl")


"""
    LAlgebra(M::Matrix)
    LAlgebra(M::Matrix ; check::Bool = false)

The L-algebra with elements {1,…,n} and multiplication table
given by M, i.e. i*j = M[i,j].

Check has a default value false. 
Set `check = true` to [`check_LAlgebra(M)`](@ref)
before constructing the L-algebra.

See also [`check_LAlgebra`](@ref), [`LAlgebraElem`](@ref).

# Examples
```jldoctest
julia> LAlgebra([2 2; 1 2])
LAlgebra([2 2; 1 2])

julia> LAlgebra([2 2; 1 2], check=true)
LAlgebra([2 2; 1 2])

julia> LAlgebra([1 1; 1 2], check=true)
ERROR: [2 2; 1 2] does not define an L-algebra
[...]
```
"""
mutable struct LAlgebra 
    matrix::Matrix
    LAlgebra(matrix::Matrix ; check::Bool = false) = check ? 
        ( check_LAlgebra(matrix) ? 
            new(matrix) : error("$M does not define an L-algebra") ) : new(matrix)
end # LAlgebra(matrix, check = true) if you want to check the matrix

"""
    LAlgebraElem(M::Matrix, value::Int)
    LAlgebra(M::Matrix ; check::Bool = false)

The element number `value` in the L-algebra defined by the matrix `M`.
See also [`LAlgebra`](@ref).

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> LAlgebraElem(A, 2)
LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
mutable struct LAlgebraElem
    algebra::LAlgebra
    value::Int
end

function Base.hash(x::LAlgebraElem, h::UInt)
    b = 0xa4e1b6fd78a06458%UInt 
    # choosen using https://www.random.org/cgi-bin/randbyte?nbytes=8&format=h
    return xor(b, xor(hash(x.value, h), h))
end


end
