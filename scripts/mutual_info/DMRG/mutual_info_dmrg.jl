
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

# function mutual_info(A,B,AB,psi)
#     Aexpt = expectation_value(psi,A)
#     Bexpt = expectation_value(psi,B)
#     ABexpt= expectation_value(psi,AB)
#     # log ( <AB>/<A>*<B>)
#     return log(ABexpt/Aexpt/Bexpt)
# end

# function Information_vec(n,psi_t,sites)
#     #Builds P_A, P_B
#     P_A = P1(1,sites)
#     P_B = P1(n+1,sites)
#     for ni in 2:n
#         P_A = apply(P_A,P1(ni,sites))
#         P_B = apply(P_B,P1(ni+n,sites))
#     end
#     #Builds P_AB
#     P_AB = apply(P_A,P_B)
#     #calculates info as a function of time
#     IAB_t = [mutual_info(P_A,P_B,P_AB,pt) for pt in psi_t]
#     return real.(IAB_t)
# end



##Constants
#Lattice Size
L =100
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

fullZi = prod(Zi)
fullPi = prod(Pi)

#time
tau = 0.01
t = [ti for ti in 0:tau:2]


 #tranverse fields strength
 h= 0
 #Hamiltonians 
 #H = Integrable_TFIM(J,h,L,sites);
H = Unintegrable_TFIM(J,J2,h,L,sites)
psi0 = random_mps(sites; linkdims=2);
energy,psi = dmrg(H,psi0;nsweeps,maxdim);

#psi = apply(fullZi,psi)
#psi /= sqrt(inner(psi',psi))
x =psi
inner(psi',apply(op("Xproj",sites[1]),psi))
x+= psi 
normalize!(x)

inner(x',apply(op("Xproj",sites[1]),x))
hi_vals = [i for i in 0:.1:2]
groundstates =[]
for hi in hi_vals
    #Hi = Integrable_TFIM(J,hi,L,sites);
    Hi = Unintegrable_TFIM(J,J2,hi,L,sites);
    energy,psi0 = dmrg(Hi,x;nsweeps,maxdim);
    push!(groundstates,psi0)
end





P1A = op("Xproj",sites[1])
P1B = op("Xproj",sites[2])
P2 = P1A*P1B
P4 = prod(Pi[1:4])
P2A= prod(Pi[1:2]); P2B = prod(Pi[3:4])

P6 = prod(Pi[1:6])
P3A = prod(Pi[1:3]); P3B = prod(Pi[4:6])
data1 =[]
data2 =[]
data3 =[]
for gs in groundstates1[2:end]

    numerator = inner(gs',apply(P2,gs))
    denom = inner(gs',apply(P1A,gs))^2
    push!(data1,log(numerator/denom))

    numerator = inner(gs',apply(P4,gs))
    denom = inner(gs',apply(P2A,gs))^2
    push!(data2,log(numerator/denom))

    numerator = inner(gs',apply(P6,gs))
    denom = inner(gs',apply(P3A,gs))^2
    push!(data3,log(numerator/denom))
end

x  = groundstates1[2]

w = inner(x',apply(P2,x))
z =inner(x',apply(P1A,x))

log(w/(z^2))

# for gs in groundstates2[2:end]
#     numerator = inner(gs',apply(P2,gs))
#     denom = inner(gs',apply(P1A,gs))^2
#     push!(data2,log(numerator/denom))
# end

plot(hi_vals[2:end],(data1),size=(600,500),dpi=300,label ="l = 1")
plot!(hi_vals[2:end],(data2),size=(600,500),dpi=300,label ="l = 2")
plot!(hi_vals[2:end],(data3),size=(600,500),dpi=300,label ="l = 3")
vline!([1])

#plot!(hi_vals[2:end],(data1),label="<Psi_0|+>=0")
#plot!(hi_vals[2:end],(data1+data2)/2,label="Average")

plot!(framestyle=:box,xlabel="h",ylabel="I(A;B)")

plot!(title="Mutual information for subsystem size")
using LaTeXStrings

inner(psi',apply(op("Xproj",sites[1]),psi))

###NORMAL TIME evolution
hi_vals = [hi for hi in 0:.1:2]
gates = [TFIM_evolution_gates(J,hi,L) for hi in hi_vals] 

times = [ti for ti in 0:0.5:20]
psi_t_h = [time_evolution(psi, g, .1, 10, 1e-8) for g in gates]


data1= []
data2= []
data3= []

#[inner(pt',apply(P2,pt))/inner(pt',apply(P1A,pt))^2 for pt in psi_t_h[2]]

for i in 1:length(hi_vals)
    gs = psi_t_h[i]


    numerator = [inner(pt',apply(P2,pt))/inner(pt',apply(P1A,pt))^2 for pt in gs]
    
    push!(data1, numerator)

    numerator = [inner(pt',apply(P4,pt))/inner(pt',apply(P2A,pt))^2 for pt in gs]
    
    push!(data2, numerator)
    numerator = [inner(pt',apply(P6,pt))/inner(pt',apply(P3A,pt))^2 for pt in gs]
    
    push!(data3, numerator)
end

i = 20
hi_vals[i]
#x= [inner(pt',apply(op("Zproj",sites[1]),pt)) for pt in psi_t_h[3]]
#plot(real.(x),dpi=300)

plot(log.(real.(data1[i])),dpi=300)
plot!(log.(real.(data2[i])))
plot!(log.(real.(data3[i])))

op('Xproj',sites[1])

########ADIABATIC

τ = 1
δt = 0.01
times = [ti for ti in 0:δt:τ]

hi = 0
hf = 0.5


time_evolve_gates = [TFIM_timedep_gates(J,hi,hf,τ,ti,δt,L) for ti in times[2:end]]
H = Unintegrable_TFIM(J,J2,hi,L,sites)
psi0 = random_mps(sites; linkdims=2);
energy,psi = dmrg(H,psi0;nsweeps,maxdim)
psi = normalize!(apply(fullZi,psi))

inner(psi',apply(P1A,psi))


psi_times = [psi]

psi_t = psi 
for g in time_evolve_gates
    psi_t = apply(g,psi_t)
    normalize!(psi_t)
    push!(psi_times,psi_t)
end


Ψ0_fid = [ inner(psi',pt) for pt in psi_times]
Ψf_fidelity = [inner(psi_f',pt) for pt in psi_times]

test = [mutual_info(P2A,P1A,psi_times[i]) for i in 1:length(times)]

plot(test,dpi=300) 


energy =[ inner(pt',apply(Hf,pt)) for pt in psi_times]


plot(times,abs.(Ψ0_fid),dpi=500, label = "< Ψ0| Ψ(t)>")
plot!(times,abs.(Ψf_fidelity),label="<Ψf|Ψ(t)>")

#Energy Plot




function mutual_info(PAB,PA,psi)
    expt_PAB = real.(inner(psi',apply(PAB,psi)))
    expt_PA  = real.(inner(psi',apply(PA,psi )))
    return log(expt_PAB/expt_PA^2)
end

h_finals = [hi for hi in 0:0.1:2]


time_evolve_gates = [TFIM_timedep_gates(J,hi,hf,τ,ti,δt,L) for ti in times[2:end]]
H = Unintegrable_TFIM(J,J2,0,L,sites)
psi0 = random_mps(sites; linkdims=2);

energy,psi0 = dmrg(H,psi0;nsweeps,maxdim)
#psi0 = normalize!(apply(fullZi,psi0))
inner(psi0,apply(P2,psi0))
energy,psi1 = dmrg(H,psi0;nsweeps,maxdim)
inner(psi1,apply(P2,psi1))

psi= normalize!(psi0+psi1)
mutual_info(P6,P3A,psi)

τ = 1
δt = .01
times = [ti for ti in 0:δt:τ]
y= time()
psi_f_hf =[]
for hf in h_finals
    println(time()-y)
    println("$hf is going now")
    psi_i = x
    psi_times = [psi]
    time_evolve_gates = [TFIM_timedep_gates(J,0,hf,τ,ti,δt,L) for ti in times[2:end]]

    for g in time_evolve_gates
        psi_i = apply(g,psi_i)
        normalize!(psi_i)
        push!(psi_times,psi_i)
    end
    push!(psi_f_hf,psi_times)
end

mut_info_1 = [(mutual_info(P2A,P1A,psi_f_hf[i][end])) for i in 1:length(h_finals)]
mut_info_2 = [(mutual_info(P4,P2A,psi_f_hf[i][end])) for i in 1:length(h_finals)]
mut_info_3 = [(mutual_info(P6,P3A,psi_f_hf[i][end])) for i in 1:length(h_finals)]
plot(h_finals,mut_info_1,ylim=(0,1),xlim=(0,2),dpi=300,label="l=1")
plot!(h_finals,mut_info_2,label="l=2")
plot!(h_finals,mut_info_3,label="l=3")
vline!([1],label="")
plot!(framestyle=:box,size=(600,500))
plot!(xlabel="Final h value",ylabel="I(A;B)",title="Adiabatic from h=0")
#savefig("Adiabatic, tau=100.png")

plot(hi_vals[2:end],(data1),size=(600,500),dpi=300,label ="l = 1")
plot!(hi_vals[2:end],(data2),size=(600,500),dpi=300,label ="l = 2")
plot!(hi_vals[2:end],(data3),size=(600,500),dpi=300,label ="l = 3")
vline!([1])
