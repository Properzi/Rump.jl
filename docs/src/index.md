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
```

### order

```@docs
in(x::LAlgebraElem, a::LAlgebra)
size(a::LAlgebra)
logical_unit
rand_elem
elements
```