push!(LOAD_PATH,"../src/")
using Documenter, Rump


makedocs(sitename="Rump Documentation", format = Documenter.HTML(prettyurls = false))
