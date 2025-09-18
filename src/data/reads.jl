
function number_of_l_algebras(n)
    if n == 1 || n == 2
        return 1
    end
    cd(@__DIR__)
    if !isdir("./LA$n")
        error("L-algebras of size $n not implemented.")
    end
    cd("./LA$n")
    if isfile("vector$n.jl")
        f = open("vector$n.jl")
        res = readlines(f)
        close(f)
        r = chop(chop(res[3]))
        vect = eval(Meta.parse(r))
        return vect[length(vect)]
        
    end
end



function small_l_algebra(n,k)
    dir=pwd()
    if n == 1
        if k != 1
            error("there is only 1 L-algebra of size 1")
        end
        return l_algebra(hcat(1))
    end
    if n == 2
        if k != 1
            error("there is only 1 L-algebra of size 2")
        end
        return l_algebra([2 2; 1 2])
    end
    cd(@__DIR__)
    if !isdir("LA$n")
        error("L-algebras of size $n not implemented")
    end
    cd("./LA$n") #folder containing database of L-alg of size n
    v = open("vector$n.jl") #file with 4 lines vn = [ // vector with in position i the number of L-algebras in the file LAn_i// cumulative vector of the previous one // ];
    vec = readlines(v)[3] #reads the cumulative vector
    close(v)
    line = chop(chop(vec)) #removing unnecessary characters
    cases = eval(Meta.parse(line)) #evaluating the string as a real vector
    i = findfirst(x -> x>=k, cases) #LA_i contains the L-algebra we want!
    prev = 0
    tot = cases[length(cases)] #total number of L-algebras of size n
    if typeof(i) == Nothing
        error("There are only $tot L-algebras of size $n")
    end
    if i==1 && !isfile("LA$n"*"_$i.jl")
        cd("./LA$n"*"_1")
        h = open("disvector$n.jl")
        disvec = readlines(h)[3]
        close(h)
        disline = chop(chop(disvec))
        discases = eval(Meta.parse(disline))
        l = findfirst(x -> x>=k, discases)
        disprev = 0
        if l != 1
            disprev = discases[l-1]
        end
        g = ""
        #only needed for n = 8
        if !(isfile("LA$n"*"_1"*"_$l.jl"))
            z = ZipFile.Reader("./LA$n"*"_1"*"_$l.zip");
            f = z.files[1]
            for j in 1:(k-disprev)+1
                g = readline(f)
            end
            r = chop(chop(g))
            M = eval(Meta.parse(r))
            cd(dir)
            return l_algebra(M)
        end
        
        
        open("LA$n"*"_1"*"_$l.jl") do f
            for j in 1:(k-disprev)+1
             g = readline(f)
            end
        end
        r = chop(chop(g))
        M = eval(Meta.parse(r))
        cd(dir)
        return l_algebra(M)
    end

    if i != 1
        prev = cases[i-1]
    end
    g = ""
    #only for the case n = 8
    if !(isfile("LA$n"*"_$i.jl"))
        z = ZipFile.Reader("./LA$n"*"_$i.zip");
        f = z.files[1]
        for j in 1:(k-prev)+1
            g = readline(f)
        end
        r = chop(chop(g))
        M = eval(Meta.parse(r))
        cd(dir)
        return l_algebra(M)
    end
    open("LA$n"*"_$i.jl") do f
        for j in 1:(k-prev)+1
          g = readline(f)
        end
    end
    r = chop(chop(g))
    M = eval(Meta.parse(r))
    cd(dir)
    return l_algebra(M)
end