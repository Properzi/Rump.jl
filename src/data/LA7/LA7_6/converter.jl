using CodecBzip2
cd(@__DIR__)

read = read("filename.g.bz2")
plain = transcode(Bzip2Decompressor, read)
write = String(plain)
open("newfilename.jl", "w") do file
    write(file, write)
end

i = 2


 for i in 1:34
    open("Ldiamond7_$i.jl","w") do io_out
        write(io_out, "diamond7_$i = [\n")
        open("diamond7_$i.jl","r") do io_in
            for res in eachline(io_in)
                r1 = replace(res, ',' => "")
                r2 = replace(r1, "]]" => "] , ")
                r3 = replace(r2, "[[" => "[")
                r4 = replace(r3, "]  [" => "; ")            
                write(io_out,r4*"\n")
            end
        write(io_out, "] ; \n")
        end
    end
 end

 for i in 1:34
    open("LA7_6_$i.jl","w") do io_out
        open("diamond7_$i.jl","r") do io_in
            for res in eachline(io_in)
                r = replace(res, "diamond7" => "LA_7_6")            
                write(io_out, r*"\n")
            end
        end
    end
 end