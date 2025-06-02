
# Rump
Rump is a Julia package for L-algebras (in future: cycle sets, rakcs, quandles, set-theoretic solutions to YBE, skew braces and all that). The package is developed by Carsten Dietzel, Lukas Gutbrunner, Silvia Properzi and Leandro Vendramin.

Created with [Pkg](https://github.com/JuliaLang/julia)

Features
--------
L-algebras
* prime elements, prime L-algebras, ideals and prime ideals
* direct and semidirect products

Databases
* L-algebras of size <8

Usage
--------
1. start Julia from the package directory (.Rump)
2. go to the package mode:
```julia-repl
    julia> ]
    (Rump) pkg>
```
3. activate the package:
```julia-repl
    (Rump) pkg> activate .
```
4.  exit the pkg mode with backspace

5. use the package:
```julia-repl
    julia> using Rump
```

The (provisorial) documentation is [here](https://htmlpreview.github.io/?https://github.com/Properzi/Rump.jl/blob/master/docs/build/index.html)


[![Build Status](https://github.com/properzi/Rump.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/properzi/Rump.jl/actions/workflows/CI.yml?query=branch%3Amaster)

[![Build Status](https://github.com/gap-packages/rig/workflows/CI/badge.svg?branch=master)](https://github.com/gap-packages/rig/actions?query=workflow%3ACI+branch%3Amaster)
[![Code Coverage](https://codecov.io/github/gap-packages/rig/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/rig)

Authors
-------
* [C. Dietzel](https://sites.google.com/view/carstendietzel/startseite)
* [L. Gutbrunner](https://www.linkedin.com/in/lukas-gutbrunner-b86aa5320/?originalSubdomain=be)
* [S. Properzi](https://properzi.github.io/)
* [L. Vendramin](https://leandrovendramin.org/)

News
----

<!---
Cite as
-------
If you have used Rump in the preparation of a paper please cite it as:
...
-->

Contributions
-------------
You are welcome to contribute with code, patches, ideas, testing and comments.


Links
-----
* [GAP - Groups, Algorithms, Programming - a System for Computational Discrete Algebra](https://www.gap-system.org/)   
* [OSCAR - Computer Algebra System](https://www.oscar-system.org/)
* [YangBaxter - GAP package](https://gap-packages.github.io/YangBaxter/)
