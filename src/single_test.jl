include("single_particle_sector.jl")
using Plots
using Statistics
using LinearAlgebra
#Model Parameters

h = 1/2;
J = -1;
L = 10;
β = 0
###########
H = H_bdg(h,L,J);
E = eigvals(H);
U = eigvecs(H)[:,1:L];
G =G_tfim(U);


###Determining Z



H = [H_bdg(h,L,J,boundary_condition = "ABC"),
    H_bdg(h,L,J,boundary_condition = "PBC")]
    

E = [eigvals(Hi)[L+1:2L] for Hi in H]
V = [eigvecs(Hi)[:,1:L] for Hi in H]
G = [G_tfim(U) for U in V]
F = [Gi[1:L,L+1:2L] for Gi in G]
G = [Gi[1:L,1:L] for Gi in G]
M = [I-2*(F[i]+G[i]) for i in 1:2]
##CHECK THIS
η = [(-1)^p*(-1)^L*(det(M[p+1])) for p in 0:1]

 

#η = [1,-1]
T = [ti for ti in .001:.01:1] #  temperature
β = 1 ./T
Z_T = []
for βi in β
    Z = []
    for Hi in 1:2
        E = eigvals(H[Hi])[L+1:2L]  # positive energies ε_μ
        E_vac =  sum(E)
        Zp = exp(βi * E_vac) * prod(1 .+ exp.(-2*βi .* E))
        Zp+= η[Hi]* exp(βi * E_vac) * prod(1 .+ exp.(-2*βi .* E))
        push!(Z, Zp)
    end
    Z_total = 0.5 * sum(Z)
    push!(Z_T,Z_total)
end

# Z_T = []
# for βi in β
#     Z = []
#     for Hi in 1:2
#         E = eigvals(H[Hi])[L+1:2L]
#         E_vac = -0.5 * sum(E)
#         Zp = exp(-βi * E_vac) * (
#                 prod(1 .+ exp.(-2βi .* E)) + η[Hi] * prod(1 .- exp.(-2βi .* E))
#               )
#         push!(Z, Zp)
#     end
#     Z_total = 0.5 * sum(Z)
#     push!(Z_T, Z_total)
# end

#plot(T,Z_T)

F_T = -log.(Z_T)/L.*T

plot(T,F_T)
#plot(β,F_T)


ϵ = []

for i in 1:length(β)
    Z = Z_T[i]
    βi = β[i]
    
end




















using LinearAlgebra
using ExponentialUtilities
Z = []
F = []
ϵ = []
for β in 0.001:0.01:2
    PF = [exp.(β.*sum(E[i])) for i in 1:2]
    T1 = [prod(1 .+ exp.(-2*β*E[i])) for i in 1:2]
    T2 = [η[i]* prod(1 .- exp.(-2*β*E[i])) for i in 1:2]
    z = sum(1/2 *PF.*(T1+T2))
    push!(Z,z)
    push!(F, -log(z)/L/β)
    T1 = [sum(E[i].*exp.(-β*E[i])) for i in 1:2]
    T2 = [η[i]* sum(E[i].* exp.(-β*E[i].+im*π)) for i in 1:2]
    num = real(T1+T2)
    push!(ϵ,sum(num)/2/z)
end

plot(Z)

#T2 = [η[i].* sum(E[i].* exp.(-β*E[i].+im*π)) for i in 1:2]
plot!(ϵ/L)



