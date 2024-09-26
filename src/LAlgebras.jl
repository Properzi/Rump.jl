"""
    check_LAlgebra(M::Matrix)

Check if the matrix M is a multiplication table of an L-algebra.

# Examples
```julia-repl
julia> check_LAlgebra([2 2; 1 2])
true

julia> check_LAlgebra([1 2; 1 2])
false
```
"""
function check_LAlgebra(M::Matrix) 
    if !isa(M, Matrix{Int})
        return false
    end
    
    if size(M, 1) != size(M, 2)
        return false
    end
    n = size(M, 2)
    lu = M[1,1] # Candidate for logical unit
    
    for i in 1:n, j in 1:n
        if (1 > M[i,j] || M[i,j] > n )
            return false
        end
    end

    for j in 1:n
        if (M[lu,j] != j || M[j,lu] != lu || M[j,j] != lu)
            return false

        end
        for i in j+1:n
            if (M[i,j] == lu && M[j,i] == lu)
                return false
            end
        end
    end
    for i in 1:n, j in i+1:n, k in 1:n
        if (M[M[i,j],M[i,k]] != M[M[j,i],M[j,k]])
            return false
        end
    end
    return true
end



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

"""
    upset(x::LAlgebraElem)

Given an element x of an L-algebra A,
returns a vector with of all the elements in A that are bigger than x.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> upset(a)
2-element Vector{LAlgebraElem}:
 LAlgebraElem(LAlgebra([2 2; 1 2]), 1)
 LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
function upset(x::LAlgebraElem)
    a = x.algebra
    return [y for y in a if y >= x]
end

function downset(x::LAlgebraElem)
    a = x.algebra
    return [y for y in a if y <= x]
end


function is_sharp(a::LAlgebra)
    for x in a, y in a 
        if x * y != x* (x * y)
            return false
        end
    end
    return true
end


function is_symmetric(a::LAlgebra)
    for x in a, y in a  
        if x * y == y && y * x != x
            return false
        end
    end
    return true
end



function is_abelian(a::LAlgebra)
    for x in a, y in a, z in a, t in a
        if (x * y) * ( z * t ) != ( x * z ) * ( y * t )
            return false
        end
    end
    return true
end


function is_linear(a::LAlgebra)
    for x in a, y in a
        if !(x <= y) || !( x >= y )
            return false
        end
    end
    return true
end


function is_discrete(a::LAlgebra)
    lu = logical_unit(a)
    for x in a, y in a
        if (x <= y) && ( x != y ) && (y != lu)
            return false
        end
    end
    return true
end

function is_semiregular(a::LAlgebra)
    for x in a, y in a, z in a
        if ((x * y) * z) * ((y * x) * z) != ((x * y) * z) * z
            return false
        end
    end
    return true
end



function is_regular(a::LAlgebra)
    if !(is_semiregular(a))
        return false
    end
    for x in a, y in a
        if (x <= y) && count(z * x == y for z in a) == 0
            return false
        end
    end
    return true
end

function is_hilbert(a::LAlgebra)
    for x in a, y in a, z in a
        if (x * (y * z)) != (x * y) * (x * z)
            return false
        end
    end
    return true
end

function is_dualBCK(a::LAlgebra)
    for x in a, y in a, z in a
        if (x * (y * z)) != (y * x) * z
            return false
        end
    end
    return true
end

function is_KL(a::LAlgebra)
    for x in a, y in a
        if !(x <= (y * x))
            return false
        end
    end
    return true
end

function is_CL(a::LAlgebra)
    lu = logical_unit(a)
    for x in a, y in a, z in a
        if (x * (y * z)) * (y * (x * z)) != lu
            return false
        end
    end
    return true
end

function is_prime_element(p::LAlgebraElem) #logical unit is not prime
    a = p.algebra
    lu = logical_unit(a)
    if p == lu
        return false
    end
    for x in a
        if !(x <= p || x * p == p)
        return false
        end
    end
    return true
end

function prime_elements(a::LAlgebra)
    E = elements(a)
    return filter(x -> is_prime_element(x), E)
end

function is_prime(a::LAlgebra)
    lu = logical_unit(a)
    for x in a
        if x!=lu && !is_prime_element(x)
            return false
        end
    end
    return true
end

function is_subLalgebra(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
    if !(issubset(s,a))
        return error("not a subset")
    end
    lu = logical_unit(a)
    if !(lu in s)
        return false
    end
    for x in s, y in s
        if !(x * y in s)
            return false
        end
    end
    return true
end


function is_invariant(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
    if !(issubset(s,a))
        return error("not a subset")
    end
    for x in a, y in s
        if !(x * y in s)
            return false
        end
    end
    return true
end

function is_ideal(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
    if !(issubset(s,a))
        return error("not a subset")
    end
    lu = logical_unit(a)
    if !(lu in s)
        return false
    end
    for x in s, y in a
        if x * y in s && !(y in s)
            return false
        end
        if !((x * y) * y in s)
            return false
        end
        if !(y * x in s || y * (x * y) in s)
            return false
        end
    end
    return true
end

function subLalgebra_generated_by(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
    T = Set([logical_unit(a)])
    S = Set(S)
    if length(S) == 0
        return T
    end
    while length(S) != 0
        union!(T,S)
        S = setdiff(union(T*S,S*T), T)
    end
    return T
end

function ideal_generated_by(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
    T = Set([logical_unit(a)])
    S = Set(S)
    L = elements(a)
    if length(S) == 0
        return T
    end
    while length(S) != 0
        union!(T,S)
        S1 = copy(S)
        for y in L
        S1 = union(S1,(S*y)*y, y*(y*S))
        end
        S = union(S1, L*S) # L*S contains S
        for x in S, y in setdiff(L, S)
            if x * y in S 
                push!(S, y)
            end
        end
        setdiff!(S, T)
    end
    return T
end

function ideal_generated_by(x::LAlgebraElem, a::LAlgebra)
    return ideal_generated_by([x],a)
end