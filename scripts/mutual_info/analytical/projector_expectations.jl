@doc """ Generates a file of data for analytical projections, so that we can compare with our other code.
"""

###Data Reading and Writing
base_folder = dirname(@__DIR__)
#Analytical Functions
include( joinpath(@__DIR__,"..","..","..","src/analytical_TFIM.jl"))


L= 8
k = [2*pi*(n+1/2)/L for n in 0:L-1];

J = 1
h1 = 0
h2 = 0
h_i = [i for i in 0:.01:5]
###
P1 = []
P2 = []
P3 = []
P4 = []
P5 = []
P6 = []
P7 = []
P8 = []

for hi in h_i
    h1 = hi
    push!(P1,P_n(1,0))
    push!(P2,P_n(2,0))
    push!(P3,P_n(3,0))
    push!(P4,P_n(4,0))
    push!(P5,P_n(5,0))
    # push!(P6,P_n(6,0))
    # push!(P7,P_n(7,0))
    # push!(P8,P_n(8,0))

end

x = [h_i,P1,P2,P3,P4]
open("test.txt","w") do io
    writedlm(io,x)
end
