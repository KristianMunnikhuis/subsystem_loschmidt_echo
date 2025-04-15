"""
Functions for using ITensor to calculate the TFIM.
"""
using ITensors, ITensorMPS
using Plots
#ITensors.compile()

######Initalization########
#Defining an onsite projector
ITensors.op(::OpName"Xproj",::SiteType"S=1/2") = 
    [1/2 1/2
     1/2 1/2]
##Constants

######Functions######
function PBC(i)
    #PBC(i+L) = i
    return mod(i-1,L)+1
end
function Integrable_TFIM(J,h)
    #Open Boundary Conditions
    os = OpSum()
    for j=1:L-1
        #sigma_x_j * sigma_x_j+1 
        os+= -J,"X",j,"X",PBC(j+1);
        os+= h, "Z", j;
    end
    os+= h, "Z", L;
    #Periodic Boundary Conditions (Manually Connect L to 1)
    os+= -J,"X",L,"X",1
    H = MPO(os, sites);
    return H
end 

function DMRG_ground_state(nsweeps, maxdim, cutoff, H)
    #Initial state is random
    psi0 = random_mps(sites; linkdims=2);
    #DMRG to find ground state
    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff);
    return psi
    # else

    #     return psi;
    # end
end

function TFIM_evolution_gates(J,h,L)
    gates = ITensor[]
    #Periodic Boundary Conditions
    for j in 1:L
        #Site 1
        s1 = sites[j]
        #Site 2
        s2 = sites[PBC(j+1)]
        #Evolve gate
        hj = -J*op("X",s1)*op("X",s2)+h*op("Z",s1)*op("I",s2)
        Gj = exp(-im*tau/2*hj)
        push!(gates,Gj)
    end
    #Append Gates
    append!(gates,reverse(gates))
    return gates
end

function time_evolution(psi, gates, tau, ttotal, cutoff)
    #initialize [Psi(t)]
    psi_times = [psi]
    #initialize psi(0)
    psi_t = psi
    for t in 0.0:tau:ttotal
        t≈ttotal && break
        #Apply time evolution operator
        psi_t = apply(gates, psi_t;cutoff)
        #normalize at each time step
        normalize!(psi_t)
        #append list
        push!(psi_times,psi_t)
    end
    #Return list
    return psi_times
end

# ####Correlators####
# function Connected_Correlator(A,B,psi)
#     expt_AB = inner(psi'*A,B*psi)
#     expt_A  = inner(psi',A,psi)
#     expt_B  = inner(psi',B,psi)
#     return 2*(expt_AB-expt_A*expt_B)
# end

# function Time_Correlator(A,B,psi_t)
#     correlations = []
#     for t in 1:length(psi_t)
#         CC = Connected_Correlator(A,B,psi_t[t])
#         push!(correlations,abs(CC))
#     end
#     return correlations
# end

##Whoole can of wormss....
# function LightCone(P,psi_t,L,starting_index,subsystem_length)
#     P_0 = P(starting_index)
#     P1_TC = []
#     indices = vcat(1:starting_index-subsystem_length,starting_index, starting_index+subsystem_length:L)
#     for i in indices
#         TC = Time_Correlator(P_0,P(i),psi_t)
#         push!(P1_TC,TC/maximum(TC)) 
        
#     end
#     plot_indices = vcat( (1:starting_index-subsystem_length).-starting_index,0, 
#                           (starting_index+subsystem_length:L).-starting_index )
#     return P1_TC, plot_indices
# end
