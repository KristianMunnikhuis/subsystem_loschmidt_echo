
include(joinpath(@__DIR__,"..","..","..","src/numerics_ITensor_TFIM.jl"))
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
L =50
J = 1;
J2 = 0;
#DMRG Parameters
cutoff_dmrg = 1e-13;
nsweeps = 20;
maxdim= [10,20,100,100,200,400,400,400,600];
data = []
#Sites
sites = siteinds("S=1/2",L)
Pi = [op("Xproj",sites[i]) for i in 1:L]
#Zi = [op("Zproj",sites[i]) for i in 1:L]
h_vals = [i for i in 0:.04:2]

function ground_state()
    
    psi0 = random_mps(sites; linkdims=2);
    while inner(psi0',apply(op("Xproj",sites[1]),psi0)) > .9
        psi0 = random_mps(sites; linkdims=2);
        H = Integrable_TFIM(J,0)
        energy,psi0 = dmrg(H,psi0;nsweeps,maxdim);
    end
    return psi0
end

psi0 = ground_state()

psi_h = [psi0]
for h in h_vals[2:end]
    H = Integrable_TFIM(J,h)
    ___,psi0 = dmrg(H, psi0; nsweeps,maxdim);
    push!(psi_h,psi0)


end

P_n_h = []
for n in 1:14
    P = prod(Pi[1:n])
    P_h = []

    for psi in psi_h
        expt_value = inner(psi',apply(P,psi))
        push!(P_h,expt_value)
    end
    push!(P_n_h,P_h)
end


plot(h_vals,P_n_h[1])
for n in 2:14
    plot!(h_vals,P_n_h[n])
end
plot!(xlims=(0,2),size=(600,500),dpi=300)