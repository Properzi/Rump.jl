# Welcome to Rump
Rump is ...

## L-algebras
```@docs
check_LAlgebra
LAlgebra
LAlgebraElem
print(x::LAlgebraElem)
*(x::LAlgebraElem, y::LAlgebraElem)
*(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}}, T::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}})
*(S::Union{Set{LAlgebraElem}, Vector{LAlgebraElem}},x::LAlgebraElem)
in(x::LAlgebraElem, a::LAlgebra)
size(a::LAlgebra)
length(a::LAlgebra)
logical_unit
rand_elem
elements
upset
downset
is_sharp
is_symmetric
is_abelian
is_linear
is_discrete
is_semiregular
is_regular
is_hilbert
is_dualBCK
is_KL
is_CL
is_prime_element
prime_elements
is_prime
is_subLalgebra
is_invariant
is_ideal
subLalgebra_generated_by
ideal_generated_by
```