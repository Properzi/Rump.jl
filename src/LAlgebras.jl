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
julia> a = LAlgebra([2 2; 1 2]);
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
    if !( x.algebra == y.algebra )
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
julia> a = LAlgebra([2 2; 1 2]);
julia> b = LAlgebra([3 1 3; 3 3 3; 1 2 3]);

julia> LAlgebraElem(a,2)*LAlgebraElem(b,2);
ERROR: elements not in the same LAlgebra

julia> LAlgebraElem(a,2)*LAlgebraElem(a,1)
LAlgebraElem(LAlgebra([2 2; 1 2]), 1)

julia> LAlgebraElem(a,1)*LAlgebraElem(a,1)*LAlgebraElem(a,1)
LAlgebraElem(LAlgebra([2 2; 1 2]), 1)

julia> LAlgebraElem(a,1)*(LAlgebraElem(a,1)*LAlgebraElem(a,1))
LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
function *(x::LAlgebraElem, y::LAlgebraElem)
    if !(x.algebra == y.algebra) 
        return error("elements not in the same L-algebra")
    end
    res = x.algebra.matrix[x.value, y.value]
    return LAlgebraElem(x.algebra, res)
end


Base.iterate(L::LAlgebra, state=1) = state > size(L) ?
     nothing : (LAlgebraElem(L,state), state+1)

"""
    in(x::LAlgebraElem, a::LAlgebra)

Determine whether `x` is an element of `a`
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
the element u in a such that x*x=x*u=u and u*x=x for all x in a.
"""
function logical_unit(a::LAlgebra)
    return LAlgebraElem(a, a.matrix[1,1])
end



function !=(x::LAlgebraElem, y::LAlgebraElem)
    if !( x.algebra == y.algebra )
        return error("elements not in the same LAlgebra")
    end
    return ( x.value != y.value)
end 

"""
    is_LAlgebraMor(domain::LAlgebra, codomain::LAlgebra, map::Vector{LAlgebraElem}) -> Bool

Check whether the given map defines a morphism of L-algebras from `domain` to `codomain`.

# Arguments
- `domain::LAlgebra`: The domain L-algebra.
- `codomain::LAlgebra`: The codomain L-algebra.
- `map::Vector{LAlgebraElem}`: A vector representing the map from `domain` to `codomain`.

# Returns
- `Bool`: `true` if `map` defines an L-algebra morphism, `false` otherwise.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([2 2; 1 2]);
julia> m = [LAlgebraElem(B,1), LAlgebraElem(B,2)];
julia> is_LAlgebraMor(A, B, m)
true
"""
function is_LAlgebraMor(domain::LAlgebra, codomain::LAlgebra, map::Vector{LAlgebraElem}) 
    m = size(domain)
    n = size(codomain)
    if length(map) != m
        return false
    end

    for i in 1:m
        if !(map[i] in codomain)
            return false
        end
    end

    for i in 1:m, j in 1:m
        x = LAlgebraElem(domain, i )
        y = LAlgebraElem(domain, j)
        if map[(x*y).value] != map[i] * map[j] 
            return false
        end
    end
    return true
end


"""
    LAlgebraMor(domain::LAlgebra, codomain::LAlgebra, map::Vector{LAlgebraElem}; check=true)

Construct an L-algebra morphism from `domain` to `codomain` given by `map`.

# Arguments
- `domain::LAlgebra`: The domain L-algebra.
- `codomain::LAlgebra`: The codomain L-algebra.
- `map::Vector{LAlgebraElem}`: A vector of elements in `codomain` representing the images of elements in `domain`.
- `check::Bool`: If `true` (default), verify that `map` defines a valid morphism.

# Returns
- `LAlgebraMor`: The L-algebra morphism object.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([2 2; 1 2]);
julia> m = [LAlgebraElem(B,1), LAlgebraElem(B,2)];
julia> LAlgebraMor(A, B, m)
LAlgebraMor(LAlgebra([2 2; 1 2]), LAlgebra([2 2; 1 2]), LAlgebraElem[LAlgebraElem(LAlgebra([2 2; 1 2]), 1), LAlgebraElem(LAlgebra([2 2; 1 2]), 2)])

julia> invalid = [LAlgebraElem(B,2), LAlgebraElem(B,1)];
julia> LAlgebraMor(A, B, invalid)
ERROR:  LAlgebraElem[LAlgebraElem(LAlgebra([2 2; 1 2]), 2), LAlgebraElem(LAlgebra([2 2; 1 2]), 1)] is not a valid L-algebra morphism from LAlgebra([2 2; 1 2]) to LAlgebra([2 2; 1 2])"""
mutable struct LAlgebraMor
    domain::LAlgebra
    codomain::LAlgebra
    map::Vector{LAlgebraElem}
    function LAlgebraMor(domain::LAlgebra, codomain::LAlgebra, map::Vector{LAlgebraElem}; check::Bool=true)
        if check
            is_LAlgebraMor(domain, codomain, map) || 
                error("$map is not a valid L-algebra morphism from $domain to $codomain")
        end
        new(domain, codomain, map)
    end
end
 

"""
    <=(x::LAlgebraElem, y::LAlgebraElem)

Check if x≤y.
"""
function <=(x::LAlgebraElem, y::LAlgebraElem)
    if !( x.algebra == y.algebra )
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
    if !( x.algebra == y.algebra )
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

Return a vector of all the elements of a.
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
return a vector of all the elements in A that are bigger than or equal to x.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> upset(a)
1-element Vector{LAlgebraElem}:
 LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
function upset(x::LAlgebraElem)
    a = x.algebra
    return [y for y in a if y >= x]
end

"""
    downset(x::LAlgebraElem)

Given an element x of an L-algebra A,
return a vector of all the elements in A that are smaller than or equal to x.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> downset(x)
2-element Vector{LAlgebraElem}:
 LAlgebraElem(LAlgebra([2 2; 1 2]), 1)
 LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
function downset(x::LAlgebraElem)
    a = x.algebra
    return [y for y in a if y <= x]
end

"""
    is_sharp(a::LAlgebra)

Check if A is a sharp L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_sharp(A)
true
```
"""
function is_sharp(a::LAlgebra)
    for x in a, y in a 
        if x * y != x* (x * y)
            return false
        end
    end
    return true
end

"""
    is_symmetric(a::LAlgebra)

Check if A is a symmetric L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_symmetric(A)
true
```
"""
function is_symmetric(a::LAlgebra)
    for x in a, y in a  
        if x * y == y && y * x != x
            return false
        end
    end
    return true
end


"""
    is_abelian(a::LAlgebra)

Check if A is an abelian L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_abelian(A)
false
```
"""
function is_abelian(a::LAlgebra)
    for x in a, y in a, z in a, t in a
        if (x * y) * ( z * t ) != ( x * z ) * ( y * t )
            return false
        end
    end
    return true
end

"""
    is_linear(a::LAlgebra)

Check if A is a linear L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_linear(A)
false
```
"""
function is_linear(a::LAlgebra)
    for x in a, y in a
        if !(x <= y) && !( x >= y )
            return false
        end
    end
    return true
end

"""
    is_discrete(a::LAlgebra)

Check if A is a discrete L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_discrete(A)
false
```
"""
function is_discrete(a::LAlgebra)
    lu = logical_unit(a)
    for x in a, y in a
        if (x <= y) && ( x != y ) && (y != lu)
            return false
        end
    end
    return true
end


"""
    is_semiregular(a::LAlgebra)

Check if A is a semiregular L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_semiregular(A)
true
```
"""
function is_semiregular(a::LAlgebra)
    for x in a, y in a, z in a
        if ((x * y) * z) * ((y * x) * z) != ((x * y) * z) * z
            return false
        end
    end
    return true
end


"""
    is_regular(a::LAlgebra)

Check if A is a regular L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_regular(A)
true
```
"""
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


"""
    is_hilbert(a::LAlgebra)

Check if A is a hilbert L-algebra
# Examples
```julia-repl
julia> A = LAlgebra([2 2; 1 2]);
julia> is_hilbert(A)
true
```
"""
function is_hilbert(a::LAlgebra)
    for x in a, y in a, z in a
        if (x * (y * z)) != (x * y) * (x * z)
            return false
        end
    end
    return true
end

"""
    is_dualBCK(a::LAlgebra)

Check if A is a dualBCK-algebra
# Examples
```julia-repl
julia> A = Lalgebra([2 2; 1 2]);
julia> is_dualBCK(A)
false
```
"""
function is_dualBCK(a::LAlgebra)
    for x in a, y in a, z in a
        if (x * (y * z)) != (y * x) * z
            return false
        end
    end
    return true
end

"""
    is_KL(a::LAlgebra)

Check if a is a KL-algebra
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> is_KL(a)
true
```
"""
function is_KL(a::LAlgebra)
    for x in a, y in a
        if !(x <= (y * x))
            return false
        end
    end
    return true
end

"""
    is_CL(a::LAlgebra)

Check if a is a CL-algebra
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> is_CL(a)
true
```
"""
function is_CL(a::LAlgebra)
    lu = logical_unit(a)
    for x in a, y in a, z in a
        if (x * (y * z)) * (y * (x * z)) != lu
            return false
        end
    end
    return true
end

"""
    is_prime_element(p::LAlgebraElem) 

Check if p is a prime element of the L-algebra.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> p = LAlgebraElem(A,1)
true
julia> n = LAlgebraElem(A,2)
false
```
"""
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

"""
    prime_elements(a::LAlgebra) 

Return a vector with all the prime elements of a 
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> prime_elements(a)
1-element Vector{LAlgebraElem}:
 LAlgebraElem(LAlgebra([2 2; 1 2]), 1)

julia> b = LAlgebra([3 1 3; 3 3 3; 1 2 3]);
julia> prime_elements(b)
1-element Vector{LAlgebraElem}:
 LAlgebraElem(LAlgebra([3 1 3; 3 3 3; 1 2 3]), 1)
```
"""
function prime_elements(a::LAlgebra)
    E = elements(a)
    return filter(x -> is_prime_element(x), E)
end

"""
    is_prime(a::LAlgebra)

Check if a is a prime L-algebra
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> is_prime(a)
true

julia> b = LAlgebra([3 1 3; 3 3 3; 1 2 3]);
julia> is_prime(b)
false
```
"""
function is_prime(a::LAlgebra)
    lu = logical_unit(a)
    for x in a
        if x!=lu && !is_prime_element(x)
            return false
        end
    end
    return true
end

"""
    is_subLalgebra(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
 
Check if s (either a subset or a vector of elements of a),
     is a L-subalgebra of a.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> is_subLalgebra([LAlgebraElem(a,1)],a)
false
julia> is_subLalgebra([LAlgebraElem(a,2)],a)
true
```
"""
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

"""
    is_invariant(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
 
Check if s (either a subset or a vector of eleemtns of a),
     is closed under the operation of a.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> is_invariant([LAlgebraElem(a,1)],a)
false
julia> is_invariant([LAlgebraElem(a,2)],a)
true
```
"""
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

"""
    is_ideal(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
 
Check if s (either a subset or a vector of eleemtns of a),
     is an ideal of a.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> is_ideal([LAlgebraElem(a,1)],a)
false
julia> is_ideal([LAlgebraElem(a,2)],a)
true
```
"""
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
        if !(y * x in s && y * (x * y) in s)
            return false
        end
    end
    return true
end

"""
    subLalgebra_generated_by(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
 
Return the L-subalgebra of a generated by the subset (or the elements in the vector) s.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> subLalgebra_generated_by([LAlgebraElem(a,1)],a)
Set{LAlgebraElem} with 2 elements:
  LAlgebraElem(LAlgebra([2 2; 1 2]), 1)
  LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
julia> subLalgebra_generated_by([LAlgebraElem(a,2)],a)
Set{LAlgebraElem} with 1 element:
  LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
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

"""
    ideal_generated_by(s::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
 
Return the ideal of a generated by the subset (or the elements in the vector) s.
# Examples
```julia-repl
julia> a = LAlgebra([2 2; 1 2]);
julia> ideal_generated_by([LAlgebraElem(a,1)],a)
Set{LAlgebraElem} with 2 elements:
  LAlgebraElem(LAlgebra([2 2; 1 2]), 1)
  LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
julia> ideal_generated_by([LAlgebraElem(a,2)],a)
Set{LAlgebraElem} with 1 element:
  LAlgebraElem(LAlgebra([2 2; 1 2]), 2)
```
"""
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
        union!(S1,(S*y)*y, y*(y*S))
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

function  ideals(a::LAlgebra)
    n = size(a)
    l = []
    for S in subsets(elements(a))
        if is_ideal(S,a)
        push!(l,S)
        end
    end
    return l
end

"""
    onlyvalues_ideals(a::LAlgebra) -> Vector{Vector{Int}}

Return all ideals of the L-algebra `a`, represented by their element values.

# Arguments
- `a::LAlgebra`: The L-algebra.

# Returns
- `Vector{Vector{Int}}`: All ideals by values.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> onlyvalues_ideals(A)
2-element Vector{Vector{Int64}}:
 [2]
 [1, 2]
"""
function onlyvalues_ideals(a::LAlgebra)
    S =  ideals(a)
    return [[x[i].value for i in 1:length(x)] for x in S]
end

"""
    ⋅(I::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}},J::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}})

Return the LAlgbera product on ideals I⋅J={a∈X : ⟨a⟩∩I ⊆ J}.
"""
function ⋅(I::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}},J::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}})
    res = Set{LAlgebraElem}()  
    if !(allequal([I[i].algebra for i in 1:length(I)]))
        return error("elements in $I are not all in the same l-algebra")
    end
    if !(allequal([J[i].algebra for i in 1:length(J)]))
        return error("elements in $J are not all in the same l-algebra")
    end
    if isempty(I)
        return error("$I is empty")
    end

    if isempty(J)
        return error("$J is empty")
    end

    a = I[1].algebra
    if !(a==J[1].algebra)
        return error("$I and $J are not in the same l-algebra")
    end

    if !is_ideal(I,a)
        return error("$I in not an ideal of $a")
    end

    if !is_ideal(J,a)
        return error("$J in not an ideal of $a")
    end

    for x in a
        if issubset((ideal_generated_by(x,a) ∩ I), J)
        push!(res, x)
        end
    end
    return res
end 

function is_prime_ideal(p::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
    if !is_ideal(p,a) || length(p) == size(a)
        return false
    end
    S =  ideals(a)
    for I in S
        if !issubset(I,p) 
            if !issubset(I ⋅ p,p)
                return false
            end
        end
    end
    return true
end

																				"""
    spec(a::LAlgebra) -> Vector{Set{LAlgebraElem}}

Return the set of all prime ideals of the L-algebra `a`.

# Arguments
- `a::LAlgebra`: The L-algebra.

# Returns
- `Vector{Set{LAlgebraElem}}`: All prime ideals.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> spec(A)
1-element Vector{Any}:
 LAlgebraElem[LAlgebraElem(LAlgebra([2 2; 1 2]), 2)]
"""
function spec(a::LAlgebra)
    filter(p->is_prime_ideal(p,a), ideals(a))
end


"""
    onlyvalues_spec(a::LAlgebra) -> Vector{Vector{Int}}

Return the set of all prime ideals of `a`, represented by their element values.

# Arguments
- `a::LAlgebra`: The L-algebra.

# Returns
- `Vector{Vector{Int}}`: All prime ideals by values.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> onlyvalues_spec(A)
1-element Vector{Vector{Int64}}:
 [2]
"""
function onlyvalues_spec(a::LAlgebra)
    S = spec(a)
    return [[x[i].value for i in 1:length(x)] for x in S]
end


"""
    dirprod(a::LAlgebra, b::LAlgebra) -> LAlgebra

Return the direct product of two L-algebras `a` and `b`.

# Arguments
- `a::LAlgebra`: The first L-algebra.
- `b::LAlgebra`: The second L-algebra.

# Returns
- `LAlgebra`: The direct product of `a` and `b`.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([2 2; 1 2]);
julia> dirprod(A,B)
LAlgebra([4 4 4 4; 3 4 3 4; 2 2 4 4; 1 2 3 4])
"""
function dirprod(a::LAlgebra,b::LAlgebra)
    n = size(a)
    m = size(b)
    s = m * n
    if n == 1 return b 
    end
    if m == 1 return a
    end
    A = zeros(Int,s,s)
    for a1 in a, a2 in a, b1 in b, b2 in b
        c = (a1*a2).value
        d = (b1*b2).value
        A[m*(a1.value -1)+b1.value ,m*(a2.value -1)+b2.value ] = m*(c-1)+d
    end
    return LAlgebra(A) #do we want it in the normal_form?
end

"""
    normalized_dirprod(a::LAlgebra, b::LAlgebra) -> LAlgebra

Return the direct product of two L-algebras `a` and `b`, then normalize it
so that the resulting L-algebra is in canonical order (normal form).

# Arguments
- `a::LAlgebra`: The first L-algebra.
- `b::LAlgebra`: The second L-algebra.

# Returns
- `LAlgebra`: The normalized direct product.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([3 3 3; 1 3 3; 1 2 3])
julia> dirprod(A,B)
LAlgebra([ 6  6  6  6  6  6; 4  6  6  4  6  6;; 4  5  6  4  5  6; 3  3  3  6  6  6; 1  3  3  4  6  6; 1  2  3  4  5  6])
julia> normalized_dirprod(A,B)
LAlgebra([ 6  6  6  6  6  6; 3  6  3  6  6  6; 4  4  6  4  6  6; 3  5  3  6  5  6; 1  4  3  4  6  6; 1  2  3  4  5  6])
"""
function normalized_dirprod(a::LAlgebra,b::LAlgebra)
    return normal_form(dirprod(a, b))
end

"""
    direct_product(a1::LAlgebra, a2::LAlgebra, ...) -> LAlgebra

Return the direct product of several L-algebras.

# Arguments
- `a1::LAlgebra, a2::LAlgebra, ...`: L-algebras to be combined.

# Returns
- `LAlgebra`: The direct product of the given L-algebras.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([2 2; 1 2]);
julia> C = LAlgebra([2 2; 1 2]);
julia> direct_product(A,B,C)
LAlgebra([8  8  8  8  8  8  8  8; 7  8  7  8  7  8  7  8; 6  6  8  8  6  6  8  8; 5  6  7  8  5  6  7  8; 4  4  4  4  8  8  8  8; 3  4  3  4  7  8  7  8; 2  2  4  4  6  6  8  8; 1  2  3  4  5  6  7  8])
"""
function direct_product((arg::LAlgebra)...)
    LA = LAlgebra([1;;])
    for (i,a) in enumerate(arg)
        LA = dirprod(LA,a)
    end
    return LA #do we want it in the normal_form?
end

"""
    normalized_direct_prod(a1::LAlgebra, a2::LAlgebra, ...) -> LAlgebra

Return the direct product of several L-algebras `a1` and `a2`..., then normalize it
so that the resulting L-algebra is in canonical order (normal form).

# Arguments
- `a1::LAlgebra, a2::LAlgebra, ...`: L-algebras to be combined.

# Returns
- `LAlgebra`: The normalized form of the direct product.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([2 2; 1 2]);
julia> C = LAlgebra([3 3 3; 1 3 3; 1 2 3])
julia> direct_product(A,B,C)
LAlgebra([12  12  12  12  12  12  12  12  12  12  12  12; 10  12  12  10  12  12  10  12  12  10  12  12; 10  11  12  10  11  12  10  11  12  10  11  12;  9   9   9  12  12  12   9   9   9  12  12  12;  7   9   9  10  12  12   7   9   9  10  12  12;  7   8   9  10  11  12   7   8   9  10  11  12;  6   6   6   6   6   6  12  12  12  12  12  12;  4   6   6   4   6   6  10  12  12  10  12  12;  4   5   6   4   5   6  10  11  12  10  11  12;  3   3   3   6   6   6   9   9   9  12  12  12;  1   3   3   4   6   6   7   9   9  10  12  12;  1   2   3   4   5   6   7   8   9  10  11  12])
julia> normalized_direct_product(A,B,C)
LAlgebra([12  12  12  12  12  12  12  12  12  12  12  12; 10  12  12  10  12  12  10  12  12  10  12  12; 10  11  12  10  11  12  10  11  12  10  11  12;  9   9   9  12  12  12   9   9   9  12  12  12;  7   9   9  10  12  12   7   9   9  10  12  12;  7   8   9  10  11  12   7   8   9  10  11  12;  6   6   6   6   6   6  12  12  12  12  12  12;  4   6   6   4   6   6  10  12  12  10  12  12;  4   5   6   4   5   6  10  11  12  10  11  12;  3   3   3   6   6   6   9   9   9  12  12  12;  1   3   3   4   6   6   7   9   9  10  12  12;  1   2   3   4   5   6   7   8   9  10  11  12])
"""
function normalized_dirprod(a::LAlgebra,b::LAlgebra)
    return normal_form(dirprod(a, b))
end

"""
    endomorphisms(a::LAlgebra) -> Vector{Vector{LAlgebraElem}}

Return the set of all endomorphisms of the L-algebra `a`.

# Arguments
- `a::LAlgebra`: The L-algebra.

# Returns
- `Vector{Vector{LAlgebraElem}}`: All endomorphisms of `a`.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> endomorphisms(A)
2-element Vector{Any}:
 LAlgebraElem[LAlgebraElem(LAlgebra([2 2; 1 2]), 1), LAlgebraElem(LAlgebra([2 2; 1 2]), 2)]
 LAlgebraElem[LAlgebraElem(LAlgebra([2 2; 1 2]), 2), LAlgebraElem(LAlgebra([2 2; 1 2]), 2)]
"""
function endomorphisms(a::LAlgebra)
	m=size(a)
	en = []
	for t in Iterators.product(ntuple(_ -> 1:m, m-1)...)
		s=[i for i in t]
		append!(s,m)
		f=map(x->LAlgebraElem(a,s[x]),1:m)
	if is_LAlgebraMor(a,a,f)
		append!(en,[f])
    end	
	end
return en
end 

"""
    onlyvalues_endomorphisms(a::LAlgebra) -> Vector{Vector{Int}}

Return the set of all endomorphisms of `a` by their element values.

# Arguments
- `a::LAlgebra`: The L-algebra.

# Returns
- `Vector{Vector{Int}}`: A vector of endomorphisms represented by their values.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> onlyvalues_endomorphisms(A)
2-element Vector{Vector{Int64}}:
 [1, 2]
 [2, 2]
"""																				
function onlyvalues_endomorphisms(a::LAlgebra)
    E = endomorphisms(a)
    return [[x[i].value for i in 1:length(x)] for x in E]
end

"""
    is_action(a::LAlgebra, b::LAlgebra, rho::Vector{Vector{LAlgebraElem}}) -> Bool

Check whether `rho` defines an action of `a` on `b`.

# Arguments
- `a::LAlgebra`: The acting L-algebra.
- `b::LAlgebra`: The L-algebra being acted on.
- `rho::Vector{Vector{LAlgebraElem}}`: The action maps.

# Returns
- `Bool`: `true` if `rho` defines a valid action, `false` otherwise.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([2 2; 1 2]);
julia> rho = [[LAlgebraElem(B,1), LAlgebraElem(B,2)], [LAlgebraElem(B,1), LAlgebraElem(B,2)]];
julia> is_action(A,B,rho)
true
"""
function is_action(a::LAlgebra, b::LAlgebra, rho::Vector{Vector{LAlgebraElem}})
	n = size(a)
	m = size(b)
    k = length(rho)
	if k != n
		return error("the length of rho is $k and the size of the active l-algebra is $n")
	end
	 
	for i in 1:n
        t = rho[i]
        if !is_LAlgebraMor(b,b,rho[i])
            return error("$rho is not an action: $t is not an l-algebra morphism")
     	end
	end

	id = map(x->LAlgebraElem(b,x), 1:m)
    	u = logical_unit(a)
	if rho[u.value] != id 
       	 return error("$rho is not an action: the logical unit does not act trivially")
    	end

    for u in a, v in a, i in 1:m
        x = rho[u.value][i]
        y = rho[v.value][i]
		if rho[(u*v).value][x.value] != rho[(v*u).value][y.value]
            return error("$rho is not an action: witnesses $u and $v")
		end
	end
return true
end

																					"""
    semidirect_prod(a::LAlgebra, b::LAlgebra, rho::Vector{Vector{LAlgebraElem}}) -> LAlgebra

Return the semidirect product of `a` and `b` with respect to the action `rho`.

# Arguments
- `a::LAlgebra`: The L-algebra providing the action.
- `b::LAlgebra`: The L-algebra being acted on.
- `rho::Vector{Vector{LAlgebraElem}}`: The action maps.

# Returns
- `LAlgebra`: The semidirect product.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> B = LAlgebra([2 2; 1 2]);
julia> rho = [[LAlgebraElem(B,1), LAlgebraElem(B,2)], [LAlgebraElem(B,1), LAlgebraElem(B,2)]];
julia> semidirect_prod(A,B,rho)
LAlgebra([4 4 4 4; 3 4 3 4; 2 2 4 4; 1 2 3 4])
"""
function semidirect_prod(a::LAlgebra, b::LAlgebra, rho::Vector{Vector{LAlgebraElem}})
	n = size(a)
	m = size(b)
	s = m * n
	if !is_action(b, a, rho)
		return error("$rho is not an action")
	end
	M = zeros(Int,s,s)
	for a1 in a, a2 in a, b1 in b, b2 in b
		l = rho[(b1*b2).value][a1.value]
		r = rho[(b2*b1).value][a2.value]
		c = (l*r).value
		d = (b1*b2).value
		M[m*(a1.value-1)+b1.value, m*(a2.value-1)+b2.value ] = m*(c-1)+d
	end
	return LAlgebra(M) #do we want it in the normal_form?
end


"""
    normal_form(a::LAlgebra) -> LAlgebra

Compute the normal form of the L-algebra `a`, where the order of indices in the resulting matrix
matches the canonical order derived from the counts of the leading unit element.

The function reorders the rows and columns of the underlying matrix of `a` so that
elements with similar structural properties (based on counts of the leading unit element)
are grouped together.

# Arguments
- `a::LAlgebra`: The input L-algebra.

# Returns
- `LAlgebra`: A new L-algebra in its normal form.

# Examples
```jldoctest
julia> A = LAlgebra([2 2; 1 2]);
julia> normal_form(A)
LAlgebra([2 2; 1 2])
julia> B = LAlgebra([2 2 3; 1 2 3; 2 2 2]);
julia> normal_form(B)
LAlgebra([3 3 3; 1 3 3; 1 2 3])
"""
function normal_form(a::LAlgebra) # order of indices fits in the order of LAlgebra
    M = a.matrix
    n = size(a)
    lu = M[1,1]

    c = [count(M[i,j] == lu for i in 1:n) for j in 1:n]
    v = sortperm(c)
    j = invperm(v)

    N = j[M[v,v]]

    return LAlgebra(N)
end
