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

#System Length
L= 10
h2 = 0

#Momentum
k = [2*pi*(n+1/2)/L for n in 0:L-1];
##Parameters
J = 1
h1 = 5
global β=1
ca(1,0)

sum([Epsilon_h(ki,h1) for ki in k])/L


dat1 = []
dat2 = []
b = range(0, stop=5, length=50)
for bi in b
    global β = bi
    push!(dat1,real(ca(0,0)))
    push!(dat2,real(ac(0,0)))
end

plot(dat1)
plot(dat2)















println(dat)
plot([i for i in 0:0.1:5],dat,xlim=(0,5),ylim=(0,.5))

#h_i range
h_i = range(0, stop=2, length=50)

b = [1/2]#range(0, stop=5, length=10)


n = [ni for ni in 1:6]
data_n = []

global β=10
BA(0,1,0)
AB(0,1,0)

for ni in n
    data = []
    for bi in b
        global β=bi
        dat = []
        for hi in h_i
            global h1 = hi

            push!(dat,P_n(ni,0))
        end
        push!(data,dat)
    end
    push!(data_n,data)
end


bi = 1
plot(h_i,abs.(data_n[1][bi]))
plot!(h_i,abs.(data_n[2][bi]))
plot!(h_i,abs.(data_n[3][bi]))
plot!(h_i,abs.(data_n[4][bi]))
plot!(h_i,abs.(data_n[5][bi]))
plot!(h_i,abs.(data_n[6][bi]))

plot!(title="$(b[bi])")
plot!(ylim=(0,.5))

# for i in 1:2:length(b) #There is some type of slight numeric error where this is equivalent to 
#     #the python code 
#    println(i)
#    plot!(h_i,real.(data[i]),label = "β = $(b[i])")
# end
# plot!(ylabel="<P>th")
# plot!(xlim=(h_i[1],h_i[end]))
x = [collect(b),collect(h_i)]
x_n = []

for data in data_n
    x_ni = copy(x)
    for d in data
        push!(x_ni,abs.(d))
    end
    #print(x_ni)
    push!(x_n,x_ni)
end

#x_n[2]-x_n[1]
#data_n[1][1]
#data_n[2][1]

for ni in n
    path = joinpath(@__DIR__, "..", "..","..", "results", "thermal_averages", "L=$(L)_n=$(ni)_TFIM_thermal_average.txt") |> normpath
    mkpath(dirname(path))

    open(path, "w") do io
        writedlm(io, x_n[ni])
    end
end
