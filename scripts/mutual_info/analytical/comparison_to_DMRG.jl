@doc """Comparison of analytical calculations of TFIM quantities
to DMRG results
"""
###Data Reading and Writing
base_folder = dirname(@__DIR__)
#Analytical Functions
include( joinpath(@__DIR__,"..","..","..","src/analytical_TFIM.jl"))




#t = [i for i in 0:0.05:T_final];

L= 14
k = [2*pi*(n+1/2)/L for n in 0:L-1];

J = 1
h1 = 0
h2 = 0
h_i = [i for i in 0:.1:2]
mut_1 = []
mut_2 = []
mut_3 = []
mut_4 = []
mut_5 = []
mut_6 = []
mut_7 = []
sub_1 = []
sub_2 = []
sub_3 = []
sub_4 = []
sub_5 = []
sub_6 = []
sub_7 = []

for hi in h_i[1:end]
    
    h1 = hi
    P_1 = P_n(1,0.0)
    P_2 = P_n(2,0.0)
    P_3 = P_n(3,0.0)
    P_4 = P_n(4,0.0)
    P_5 = P_n(5,0.0)
    P_6 = P_n(6,0.0)
    P_7 = P_n(7,0.0)
    push!(mut_1,log(P_2/P_1^2))
    push!(mut_2,log(P_4/P_2^2))
    push!(mut_3,log(P_6/P_3^2))
    push!(mut_4,log(P_n(8,0.0)/P_4^2))
    push!(mut_5,log(P_n(10,0.0)/P_5^2))
    push!(mut_6,log(P_n(12,0.0)/P_6^2))
    push!(mut_7,log(P_n(14,0.0)/P_7^2))

    push!(sub_1,P_1)
    push!(sub_2,P_2)
    push!(sub_3,P_3)
    push!(sub_4,P_4)
    push!(sub_5,P_5)
    push!(sub_6,P_6)
    push!(sub_7,P_7)
end




plot(h_i,abs.(mut_1),dpi=300, size = (600,500), title = "Mutual Information", xlabel = "h", ylabel = "I")
plot!(h_i,abs.(mut_2))
plot!(h_i,abs.(mut_3))
plot!(h_i,abs.(mut_4))
plot!(h_i,abs.(mut_5))
plot!(h_i,abs.(mut_6))
plot!(h_i,abs.(mut_7))


##################DMRG
#DMRG Functions

include( joinpath(@__DIR__,"..","..","..","src/numerics_ITensor_TFIM.jl"))

#Let's now define some Functions
function expectation_value(psi,O)
    #<psi|O|psi>
    return inner(psi',apply(O,psi))
end

##Constants
#Lattice Size
L =14
J = 1;
J2 = 0;
#DMRG Parameters
cutoff_dmrg = 1e-13;
nsweeps = 20;
maxdim= [10,20,100,100,200,400,400,400,600];
data = []
#Sites
sites = siteinds("S=1/2",L);
Xi = [op("Xproj",sites[i]) for i in 1:L]
Zi = [op("Zproj",sites[i]) for i in 1:L]

fullZi = prod(Zi)

#Initialize First state
h_dmrg = [i for i in 0:.1:2]

h0 = 0

H = Integrable_TFIM(J,h0)
psir = random_mps(sites; linkdims=2);
energy, psi0 = dmrg(H,psir;nsweeps,maxdim);
psi0 = apply(fullZi,psi0)
normalize!(psi0)
#energy, psi0 = dmrg(H,psi0;nsweeps,maxdim);

psi_hi = [psi0]
for hi in h_dmrg[2:end]
    Hi = Integrable_TFIM(J,hi)
    energy,psi = dmrg(Hi,psi_hi[end];nsweeps,maxdim);
    push!(psi_hi,psi)
end

mut_1_dmrg = []
mut_2_dmrg = []
mut_3_dmrg = []
mut_4_dmrg = []
mut_5_dmrg = []
mut_6_dmrg = []
mut_7_dmrg = []
P1 = Xi[1]
P2 = prod(Xi[1:2])
P3 = prod(Xi[1:3])
P4 = prod(Xi[1:4])
P5 = prod(Xi[1:5])
P6 = prod(Xi[1:6])
P7 = prod(Xi[1:7])
P8 = prod(Xi[1:8])
P9 = prod(Xi[1:9])
P10 = prod(Xi[1:10])
P11 = prod(Xi[1:11])
P12 = prod(Xi[1:12])
P13 = prod(Xi[1:13])
P14 = prod(Xi[1:14])


for i in 1:length(psi_hi)
    pt = psi_hi[i]
    push!(mut_1_dmrg,log( expectation_value(pt,P2)/expectation_value(pt,P1)^2   ))
    push!(mut_2_dmrg,log( expectation_value(pt,P4)/expectation_value(pt,P2)^2   ))
    push!(mut_3_dmrg,log( expectation_value(pt,P6)/expectation_value(pt,P3)^2   ))
    push!(mut_4_dmrg,log( expectation_value(pt,P8)/expectation_value(pt,P4)^2   ))
    push!(mut_5_dmrg,log( expectation_value(pt,P10)/expectation_value(pt,P5)^2   ))
    push!(mut_6_dmrg,log( expectation_value(pt,P12)/expectation_value(pt,P6)^2   ))
    push!(mut_7_dmrg,log( expectation_value(pt,P14)/expectation_value(pt,P7)^2   ))
end


scatter(h_dmrg[2:end],mut_1_dmrg[2:end],ms= 5, label = "DMRG l = 1",color="blue")
scatter!(h_dmrg[2:end],mut_2_dmrg[2:end],ms= 5, label = "DMRG l = 2", color = "red")

scatter!(h_dmrg[2:end],mut_3_dmrg[2:end],ms= 5, label = "DMRG l = 3", color = "green")
scatter!(h_dmrg[2:end],mut_4_dmrg[2:end],ms= 5, label = "DMRG l = 4", color = "purple")
scatter!(h_dmrg[2:end],mut_5_dmrg[2:end],ms= 5, label = "DMRG l = 4", color = "aquamarine")
scatter!(h_dmrg[2:end],mut_6_dmrg[2:end],ms= 5, label = "DMRG l = 4", color = "orchid")
scatter!(h_dmrg[2:end],mut_7_dmrg[2:end],ms= 5, label = "DMRG l = 4", color = "orange")

plot!(h_i,real.(mut_1), lw = 3,color = "blue",label = "l = 1" )
plot!(h_i,real.(mut_2), lw = 3, color = "red", label = "l = 2")
plot!(h_i,real.(mut_3), lw = 3, color = "green", label = "l = 3")
plot!(h_i,real.(mut_4), lw = 3, color = "purple", label = "l = 4")
plot!(h_i,real.(mut_5), lw = 3, color = "aquamarine", label = "l = 5")
plot!(h_i,real.(mut_6), lw = 3, color = "orchid", label = "l = 6")
plot!(h_i,real.(mut_7), lw = 3, color = "orange", label = "l = 7")

plot!(xlim=(0,2),ylim=(0,1),framestyle=:box, size =(600,500),legend= :bottomleft)
plot!(xlabel = "h", ylabel = "I(A;B)")

savefig("DMRG_Analytical_Comparison.png")

function derivative(f, x)
    N = length(x)
    df = similar(f)
    df[1] = (f[2] - f[1]) / (x[2] - x[1])                    # forward diff
    for i in 2:N-1
        df[i] = (f[i+1] - f[i-1]) / (x[i+1] - x[i-1])        # central diff
    end
    df[N] = (f[N] - f[N-1]) / (x[N] - x[N-1])                # backward diff
    return df
end


# function derivative(f,x)
#     derivative = [(f[i+1]-f[i])/(x[i+1]-x[i]) for i in 1:length(x)-1]
#     return derivative
# end

plot(h_i,derivative(real.(mut_1),h_i))
scatter!(h_dmrg,derivative(real.(mut_1_dmrg),h_dmrg))