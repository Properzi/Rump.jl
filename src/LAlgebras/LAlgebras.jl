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
ERROR: [1 1; 1 2] does not define an L-algebra
[...]
```
"""
mutable struct LAlgebra 
    matrix::Matrix
    LAlgebra(matrix::Matrix ; check::Bool = false) = check ? 
        ( check_LAlgebra(matrix) ? 
            new(matrix) : error("$matrix does not define an L-algebra") ) : new(matrix)
end # LAlgebra(matrix, check = true) if you want to check the matrix

function ==(a::LAlgebra, b::LAlgebra) # when?
     return (a.matrix == b.matrix )
end

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


"""
    print(x::LAlgebraElem)

When applied to L-algebra elements,
 prints (using println) the  value of x.
"""
function print(x::LAlgebraElem)
    println(x.value)
end




function ==(x::LAlgebraElem, y::LAlgebraElem)
    if ( x.algebra !== y.algebra )
        # return false
        return error("elements not in the same L-algebra")
    end
    return ( x.value == y.value)
end

"""
    *(x::LAlgebraElem, y::LAlgebraElem)

When applied to L-algebra elements of the same L-algebra,
multiplies them in the L-algebra.
# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);

julia> LAlgebraElem(A,2)*LAlgebraElem(A,1)
LAlgebraElem(LAlgebra([2 2; 1 2]), 1)

julia> LAlgebraElem(A,1)*LAlgebraElem(A,1)*LAlgebraElem(A,1)
LAlgebraElem(LAlgebra([2 2; 1 2]), 1)

julia> LAlgebraElem(A,1)*(LAlgebraElem(A,1)*LAlgebraElem(A,1))
LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
function *(x::LAlgebraElem, y::LAlgebraElem)
    if x.algebra != y.algebra 
        return error("elements not in the same LAlgebra")
    end
    res = x.algebra.matrix[x.value, y.value]
    return LAlgebraElem(x.algebra, res)
end


Base.iterate(L::LAlgebra, state=1) = state > size(L) ?
     nothing : (LAlgebraElem(L,state), state+1)

"""
    in(x::LAlgebraElem, a::LAlgebra)

Determine wether `x` is an element of `a`
in the sense that `x.algebra==a`
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([3 1 3; 1 3 3; 1 2 3]);
julia> LAlgebraElem(A,2) in A
true

julia> LAlgebraElem(B,2) in A
false
```
"""
function in(x::LAlgebraElem, a::LAlgebra)
    return (x.algebra == a)
end

"""
    size(a::LAlgebra)

When applied to an L-algebra returns its cardinality.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> size(A)
2
```
"""
function size(a::LAlgebra)
    return size(a.matrix, 2)
end


"""
    length(a::LAlgebra)
 
 When applied to an L-algebra returns its cardinality.
 
 # Examples
 ```jldoctest
 julia> A = LAlgebra([2 2; 1 2]);
 julia> length(A)
 2
 ```
 """ 
Base.length(L::LAlgebra) = size(L)

"""
    logical_unit(a::LAlgebra))

Return the logical unit of the L-algebra a, i.e.
the element u in a such that x*x=x*u=1 and 1*x=x for all x in a.
"""
function logical_unit(a::LAlgebra)
    return LAlgebraElem(a, a.matrix[1,1])
end



function !=(x::LAlgebraElem, y::LAlgebraElem)
    if ( x.algebra !== y.algebra )
        return error("elements not in the same LAlgebra")
    end
    return ( x.value != y.value)
end 

"""
    <=(x::LAlgebraElem, y::LAlgebraElem)

Check if x≤y.
"""
function <=(x::LAlgebraElem, y::LAlgebraElem)
    if ( x.algebra !== y.algebra )
        return error("elements not in the same LAlgebra")
    end
    a = x.algebra
    lu = logical_unit(a)
    return ( x * y == lu)
end

"""
    <(x::LAlgebraElem, y::LAlgebraElem)

Check if x<y.
"""
function <(x::LAlgebraElem, y::LAlgebraElem)
    return (x <= y && x != y)
end 

"""
    >=(x::LAlgebraElem, y::LAlgebraElem)

Check if x≥y.
"""
function >=(x::LAlgebraElem, y::LAlgebraElem)
    if ( x.algebra !== y.algebra )
        return error("elements not in the same LAlgebra")
    end
    a = x.algebra
    lu = logical_unit(a)
    return ( y * x == lu)
end

"""
    >(x::LAlgebraElem, y::LAlgebraElem)

Check if x>y.
"""
function >(x::LAlgebraElem, y::LAlgebraElem)
    return (x >= y) && (x != y)
end

"""
    *(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}},x::LAlgebraElem)

Return the set of products {s*t:s∈S,t∈T}.
"""
function *(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, T::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}})
    res = Set{LAlgebraElem}()  
    for s in S, t in T
        push!(res,s*t)
    end
    return res
end

"""
    *(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}},x::LAlgebraElem)

Return the set of products {s*x∣s∈S}.
"""
function *(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}},x::LAlgebraElem)
    return S * [x]
end

"""
    *(x::LAlgebraElem,S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}})

Return the set of products {x*s∣s∈S}.
"""
function *(x::LAlgebraElem,S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}})
    return [x] * S
end

"""
    rand_elem(a::LAlgebra)

Return a random element of a.
"""
function rand_elem(a::LAlgebra)
    return LAlgebraElem(a,rand(1:size(a)))    
end

"""
    elements(a::LAlgebra)

Return a vector with of all the elements of a.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> elements(a)
2-element Vector{LAlgebraElem}:
 LAlgebraElem(LAlgebra([2 2; 1 2]), 1)
 LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
function elements(a::LAlgebra)
    n = size(a)
    return [LAlgebraElem(a,y) for y in 1:n]
end
