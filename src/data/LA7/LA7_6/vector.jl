cd(@__DIR__)
v = zeros(Int, 34)
w = zeros(Int, 34)
for i in 1:34
    v[i] =countlines("LA7_6_$i.jl")-2
    w[i] = sum(v)
end

open("diamond_vector7.jl", "w") do io
    write(io, "diamond7 = [\n $v , \n $w , \n ] ;")
end