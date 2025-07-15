
using Pkg
Pkg.add("ITensors")
Pkg.add("ITensorMPS")
Pkg.add("Plots")

include(joinpath(@__DIR__,"..","..","src/numerics_ITensor_TFIM.jl"))
##Plotting

ENV["GKS_ENCODING"] = "utf-8"
using Plots
using DelimitedFiles


#Let's now define some Functions
function expectation_value(psi,O)
    #<psi|O|psi>
    return inner(psi',O,psi)
end

##Constants
#Lattice Size
L =100
J = 1;
J2 = 0;
g  = 1;
#DMRG Parameters
cutoff_dmrg = 1e-10;
nsweeps = 20;
maxdim= [10,20,100,100,200,400,400,400,600];
data = []
#Sites
sites = siteinds("S=1/2",L)

H = Integrable_TFIM(J,g)

#Finding Ground  state 
psi0 = random_mps(sites; linkdims=2)
energy,psi0 = dmrg(H,psi0;nsweeps,maxdim)

x=correlation_matrix(psi0,"X","X")

plot(x[1,:])

f= psi0
n = [i for i in 1:50]
N =inner(f', apply(op("X",sites[1])*op("X",sites[2]),f))
dat = []
for i in n[4:50]

    push!(dat,inner(f',apply(op("X",sites[1])*op("X",sites[2])*op("X",sites[i])*op("X",sites[i+1]),f))- N^2)
end

plot(n[4:50],dat)
plot!(n,n.^-2)
plot!(xscale=:log10,yscale=:log10)
xlims!(4,25)
n = [i for i in 1:50]
plot(n,x[1,2:51])
plot!(n,x[1,2].*n.^(-1/4))
plot!(xscale=:log10)
plot!(yscale=:log10)
xlims!(1,50)
##SCALING OF XX 
#Let's now define some Functions
function expectation_value(psi,O)
    #<psi|O|psi>
    return inner(psi',O,psi)
end




for i in n[2:50]
    X4i = op("X",sites[1])*op("X",sites[2])*op("X",sites[i])*op("X",sites[i+1])
    X2 = op("X",sites[1])*op("X",sites[2])

    value = inner(psi0',apply(X4i,psi0))-inner(psi0',apply(X2,psi0))^(2)
    push!(dat,value)
end

plot(n[2:50],abs.(dat))