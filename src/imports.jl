
#OSCAR packages?
#GAP enviroment
using IterTools

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
