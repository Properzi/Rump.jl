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
    if x.algebra != y.algebra 
        return error("elements not in the same L-algebra")
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
    is_LAlgebraMor(M::Matrix)

Check if 
# Examples
```julia-repl
julia> 
true

julia> 
false
```
"""
function is_LAlgebraMor(domain::LAlgebra, codomain::LAlgebra, map::Vector{LAlgebraElem}) 
    m = size(domain)
    n = size(codomain)
    if length(f) != m
        return false
    end

    for i in 1:m
        if !(f[i] in codomain)
            return false
        end
    end

    for i in 1:m, j in 1:m
        x = LAlgebraElem(domain, i )
        y = LAlgebraElem(domain, j)
        if f[(x*y).value] != f[i] * f[j] 
            return false
        end
    end
    return true
end


"""
    LAlgebraMor(a::LAlgebra, b::LAlgebra, map::Vector{LAlgebraElem})
    LAlgebraMor(a::LAlgebra, b::LAlgebra, map::Vector{LAlgebraElem}, check=false)



# Examples
```jldoctest
julia> 
```
"""
mutable struct LAlgebraMor
    domain::LAlgebra
    codomain::LAlgebra
    map::Vector{LAlgebraElem}
    LAlgebra(matrix::Matrix ; check::Bool = false) = check ? 
        ( is_LAlgebraMor(domain, codomain, map) ? 
            new(map) : error("$map is not a L-algebra morphism") ) : new(map)
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
return a vector with of all the elements in A that are bigger than x.
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
return a vector with of all the elements in A that are bigger than x.
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
        if !(x <= y) || !( x >= y )
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
julia> is_dualBCK(a)
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

Check if a is a KL-algebra
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
 
Check if s (either a subset or a vector of eleemtns of a),
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
        if !(y * x in s || y * (x * y) in s)
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
 
Return the L-subalgebra of a generated by the subset (or the elements in the vector) s.
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

function spec(a::LAlgebra)
    n = size(a)
    l = []
    for S in subsets(elements(a))
        if is_ideal(S,a)
        push!(l,S)
        end
    end
    return l
end

function is_prime_ideal(p::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, a::LAlgebra)
    if !is_ideal(p,a) || length(p) == size(a)
        return false
    end
    S = spec(a)
    for I in S
        if !issubset(I,p) & !issubset(p*I,p)
            return false
        end
    end
    return true
end

function prime_spec(a::LAlgebra)
    filter(p->is_prime_ideal(p,a), spec(a))
end


function dirprod(L::LAlgebra,M::LAlgebra)
    n = size(L)
    m = size(M)
    s = m*n
    if n == 1 return M 
    end
    if m == 1 return L
    end
    A = zeros(Int,s,s)
    for a1 in L, a2 in L, b1 in M, b2 in M
        a = (a1*a2).value
        b = (b1*b2).value
        A[m*(a1.value -1)+b1.value ,m*(a2.value -1)+b2.value ] = m*(a-1)+b
    end
    return LAlgebra(A)
end


function direct_product((arg::LAlgebra)...)
    LA = LAlgebra([1;;])
    for (i,a) in enumerate(arg)
        LA = dirprod(LA,a)
    end
    return LA
end




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