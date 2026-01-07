@doc """TO BE USED WITH CLUSTER. Generates a file of data for analytical projections, so that we can compare with our other code.
"""

###Data Reading and Writing
base_folder = dirname(@__DIR__)
#Analytical Functions
include( joinpath(@__DIR__,"..","..","src/analytical_TFIM.jl"))

#Subsystem Length
L = parse(Int, ARGS[1])


#Momentum
k = [2*pi*(n+1/2)/L for n in 0:L-1];
##Parameters
J = 1
h1 = 0
h2 = 0
#h_i range
h_i = [i for i in 0:.1:5]
###
Projectors = []
for li in 1:L
    P =[]
    for hi in h_i
        global h1 = hi
        push!(P,real(P_n(li,0)))
    end
    push!(Projectors,P)

end

x = [h_i]
for P in Projectors
    println(P)
    push!(x,(P))
end

path = joinpath(@__DIR__, "..","..","results", "partial_projectors", "L=$(L)_TFIM_projectors.txt") |> normpath
mkpath(dirname(path))

open(path, "w") do io
    writedlm(io, x)
end
