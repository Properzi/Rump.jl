# Rump
[![Build Status](https://github.com/properzi/Rump.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/properzi/Rump.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://htmlpreview.github.io/?https://github.com/Properzi/Rump.jl/blob/master/docs/build/index.html)


Rump is a designed for the study of L-algebras. It provides tools for constructing L-algebras, checking their properties (e.g., primality, ideals, subalgebras), and performing algebraic operations such as direct products and semidirect products.

## Installation

Option 1: Using Pkg.add
1. Start Julia.
2. In the Julia REPL, add the package using Pkg.add:
   ```julia
   julia> using Pkg
   julia> Pkg.add(PackageSpec(url="https://github.com/Properzi/Rump.jl"))
   ```

Option 2: Manual Setup
1. Clone the repository and change to the package directory (.Rump).

## Usage

1. Start Julia from the package directory.
2. Enter the package mode by pressing `]`:
   ```julia
   julia> ]
   (Rump) pkg>
   ```

3. Activate the package environment:
   ```julia
   (Rump) pkg> activate .
   ```

4. Instantiate the package to make sure all dependencies are installed:
   ```julia
   (Rump) pkg> instantiate
   ```

5. Exit the package mode by pressing backspace.

6. Use the package:
   ```julia
   julia> using Rump
   ```

You can now use the Rump package in your Julia scripts or REPL session.




