@doc """ Generates a file of data for analytical projections, so that we can compare with our other code.
"""

###Data Reading and Writing
base_folder = dirname(@__DIR__)
#Analytical Functions
include( joinpath(@__DIR__,"../../../src/analytical_TFIM.jl"))

function nthermal(k,h)
    ϵ= Epsilon_h(k,h)
    return 1/(1+exp(β*ϵ))
end
#Define Thermal expectations
function aa(l,t)
    #<c_i c_i+l>(t)
    sum = 0
    for ki in k
        θ = Theta(ki,h1)
        ϕ = Theta(-ki,h1)
        n = nthermal(ki,h1)
        phase = exp(-im*l*ki)
        term_1= im*cos(θ/2)*sin(ϕ/2)*(1-n)
        term_2= im*sin(θ/2)*cos(ϕ/2)*n
        sum += (term_1+term_2)*phase
    end
    return sum/L
end 
function cc(l,t)
    #<c_i c_i+l>(t)
    sum = 0
    for ki in k
        θ = Theta(ki,h1)
        ϕ = Theta(-ki,h1)
        n = nthermal(ki,h1)
        ###DO THEY HAVE THE SAME PHASE?!?!?
        phase = exp(im*l*ki)
        term_1= -im*cos(θ/2)*sin(ϕ/2)*(n)
        term_2= -im*sin(θ/2)*cos(ϕ/2)*(1-n)
        sum += (term_1+term_2)*phase
    end
    return sum/L
end 

function ca(l,t)
    #<c_i c_i+l>(t)
    sum = 0
    for ki in k
        θ = Theta(ki,h1)
        ϕ = Theta(-ki,h1)
        n = nthermal(ki,h1)
        phase = exp(im*l*ki)
        term_1= cos(θ/2)^2*n
        term_2= sin(ϕ/2)^2*(1-n)
        sum += (term_1+term_2)*phase
    end
    return sum/L
end 
function ac(l,t)
    #<c_i c_i+l>(t)
    sum = 0
    for ki in k
        θ = Theta(ki,h1)
        ϕ = Theta(-ki,h1)
        n = nthermal(ki,h1)
        phase = exp(-im*l*ki)
        term_1= cos(θ/2)^2*(1-n)
        term_2= sin(ϕ/2)^2*(n)
        sum += (term_1+term_2)*phase
    end
    return sum/L
end 

#Subsystem Length
L= 8

#Momentum
k = [2*pi*(n+1/2)/L for n in 0:L-1];
##Parameters
J = 1
h1 = 0
h2 = 0
#h_i range
h_i = [i for i in 0:.01:5]
b = exp10.(range(-2,log10(3),15))


data = []
for bi in b
    global β=bi
    dat = []
    for hi in h_i
        global h1 = hi

        push!(dat,P_n(2,0))
    end
    push!(data,dat)
end
data
plot(xlabel="h")

for i in 1:2:length(b) #There is some type of slight numeric error where this is equivalent to 
    #the python code 
   println(i)
   plot!(h_i,real.(data[i]),label = "β = $(b[i])")
end
plot!(ylabel="<P>th")
plot!(xlim=(h_i[1],h_i[end]))
x=vcat(h_i,data)

x = [b,h_i]
for d in data
    push!(x,real.(d))
end
x
path = joinpath(@__DIR__, "..", "..","..", "results", "thermal_averages", "L=$(L)_TFIM_thermal_average.txt") |> normpath
mkpath(dirname(path))

open(path, "w") do io
    writedlm(io, x)
end