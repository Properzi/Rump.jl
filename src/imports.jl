#standard packages
using Pkg
using Random
using RandomExtensions
using UUIDs
using TOML

#OSCAR packages?


# import stuff from Base for which we want to provide extra methods
import Base: 
    print,
    ==,
    !=, 
    <=, 
    <, 
    >=, 
    >, 
    in, 
    size, 
    *
