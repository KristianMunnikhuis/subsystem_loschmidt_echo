using Plots
using SkewLinearAlgebra
using LinearAlgebra
using Combinatorics
##Necessary Functions to calculate everything. See my Latex Doc / the TFIM Quench Paper cited in it.


function Epsilon_0(k, h)
    #Diagonal Dispersion 
        return 2*h - 2*J*cos(k)
    end
    
    function Epsilon_h(k,h)
    #Full Dispersion
        return sqrt( Epsilon_0(k,h)^2+ 4*J*J*sin(k)^2 )
    end
    
    function Theta(k, h)
    #Diagonalization Angle
        x =-im * log( (Epsilon_0(k,h)+ im*2*J*sin(k))/Epsilon_h(k,h)  )
        return x.re
    end
    function Delta(k,h1,h2)
    #Pre and Post Quench Angle difference
        return Theta(k,h2)-Theta(k,h1)
    end
    
    function A(k, t, h1, h2)
    #Bogoliouv Amplitude
        th = Theta(k, h2)
        del = Delta(k,h1,h2)
        eps = Epsilon_h(k, h2)
        
        return cos(th / 2) * cos(del / 2) * exp(-im * eps * t) +
               sin(th / 2) * sin(del / 2) * exp( im * eps * t)
    end
    
    function B(k, t, h1, h2)
    # |A|^2 + |B|^2 = 1, but also this is its definition.
        th = Theta(k, h2)
        del = Delta(k,h1,h2)
        eps = Epsilon_h(k, h2)
        
        return im*sin(th / 2) * cos(del / 2) * exp(im * eps * t) +
               -im*cos(th / 2) * sin(del / 2) * exp(- im * eps * t)
    end
    
    ##Defining Types of two point operators. 
    
    
    ####I believe these to be fixed now. In that if you calculate out <sigma_j sigma_k> and compute linearly left to right, you can verbatum put these in 
    # and get the correct answer. Eg. 
    # sigma_x_j sigma_x_j+1 = (c_d_j - c_j) (c_d_j+1 + c_j+1) = cc(1,t) + ca(1,t) - ac(1,t) - aa(1,t) 
    ####
    # <c_i c_i+l>(t)
    function aa( l,t)
        sum = 0
        for ki in k
            Ap = A(ki,t,h1,h2)
            Bp = B(-ki,t,h1,h2)
            phase = exp(-im*l*ki)
            sum += Ap*Bp*phase
        end
        return sum/L
    end 
    #<c^dag_i c^dag_i+l> (t) 
    function cc(l,t)
        sum = 0
        for ki in k
            Ap = A(ki,t,h1,h2)
            Bp = B(-ki,t,h1,h2)
            phase = exp(-im*l*ki)
    
            sum += conj(Ap)*conj(Bp)*phase
        end
        return sum/L
    end
    #<c^dag_i c_i > (t)
    function ca(l,t)
        sum = 0
        for ki in k
            Bp = B(ki,t,h1,h2)
            B2 =  Bp*conj(Bp)
            phase = exp(im*l*ki)
            sum+= B2*phase
        end
        return sum/L
    end
    #<c_i c^dag_i> (t)
    function ac(l,t)
        sum = 0
        for ki in k
            Ap = A(ki,t,h1,h2)
            A2 =  Ap*conj(Ap)
            phase = exp(-im*l*ki)
            sum+= A2*phase
        end
        return sum/L
      end
    ###MAJORNA FERMION REPRESENTATION 
    #Since we have calcuated all expectation values of creation and annihliation operators, we can now calculate 
    #Expectaiton values of majorna fermions A and B 
    function AA(i,j,t)
        l = j-i
        return cc(l,t)+ca(l,t)+ac(l,t)+aa(l,t)
    end
    function BB(i,j,t)
        l = j-i
        return -cc(l,t)-aa(l,t)+ca(l,t)+ac(l,t)
    end
    function AB(i,j,t)
        l = j-i
        return im*(-aa(l,t)+ac(l,t)-ca(l,t)+cc(l,t)) 
    end
    function BA(i,j,t)
        l = j-i
        return im*(-aa(l,t)-ac(l,t)+ca(l,t)+cc(l,t))
    end
    
    
    # Function to construct the antisymmetric M matrix
    function construct_M_matrix(fermion_operators,  t)
        """
        fermion_operators: List of fermionic operators (site index, operator type)
                          Example: [(i, 'A'), (j, 'A'), (k, 'B'), (z, 'B')]
        correlation_function: A function that returns two-point correlators ⟨f_i f_j⟩
        """
    
        l = length(fermion_operators)  # Number of operators
        M = zeros(ComplexF64, l, l)   # Initialize antisymmetric matrix
    
        for μ in 1:l, ν in (μ+1):l
            op_1 = fermion_operators[μ]
            op_2 = fermion_operators[ν]
            M[μ, ν] = correlation_function(op_1,op_2,t)
    
    
            # Enforce antisymmetry: M_νμ = -M_μν
            M[ν, μ] = -M[μ, ν]
        end
    
        return M
    end
    
    # Example correlation function (replace with actual post-quench correlations)
    function correlation_function(op_string_1, op_string_2, t)
        #Unpack
        i, op1 = op_string_1
        j, op2 = op_string_2
        ##
        if op1 == "A" && op2 == "A"
            return AA(i,j, t)  # ⟨A_i A_j⟩
        elseif op1 == "B" && op2 == "B"
            return BB(i,j,t)  # ⟨B_i,B_j⟩
        elseif op1 == "A" && op2 == "B"
            return AB(i,j,t)  # ⟨A_i A_j⟩
        elseif op1 == "B" && op2 == "A"
            return BA(i,j,t) # ⟨B_i A_j⟩
        end
    end
    
    function sigma_general(indices,t)
        ODD = 0
        indices = sort(indices)
        function remove_duplicates_in_pairs(vec)
            unique_vals = unique(vec)
            freqs = [count(==(val), vec) for val in unique_vals]
            filtered_vals = unique_vals[freqs .% 2 .!= 0]  # Keep elements that don't fully pair up
            return filtered_vals
        end
        indices = remove_duplicates_in_pairs(indices)
        #println(indices)
        if isempty(indices)
            return 1
        end
        #If there are an odd number of sigmas then we must use cluster decomposition to approximate the answer
        if isodd(length(indices))
            ODD = 1
            l = L//2
            decomposition_indices = indices.+ l
            indices = vcat(indices,decomposition_indices)
        end
        #Build String List
        string_list = []
        for i in 1:length(indices)
            if isodd(i)
                push!(string_list,(indices[i],"B"))
            else
                push!(string_list,(indices[i],"A"))
            end
        end 
        #Inter-Term Strings
        #Inter-Term Strings are formed by pairs of sigma_x
        #Eg. <sigma_i sigma_j> has a string between i+1:j-1 
        #Since we sorted already, we can go forward in steps of 2
        pf = 1
        for i in 1:2:length(indices)
            for a in indices[i]+1:indices[i+1]-1
                pf *= -1*(im)^(indices[i+1]-indices[i])
                push!(string_list, (a,"B"))
                push!(string_list, (a,"A"))
            end
        end
        M = construct_M_matrix(string_list,t)
        if ODD == 1
            return abs(sqrt(pfaffian(M)))
        else
            return pf*pfaffian(M)
        end
    end      
            
    function binomial_expansion(indices)
        n = length(indices)
        pairs = collect(combinations(indices, 2))
        return pairs
    end
    
    function unique_elements_and_frequencies(vec)
        unique_vals = unique(vec)
        freqs = [count(==(val), vec) for val in unique_vals]
        return unique_vals, freqs
    end
    
    function P_n(n,t)
        #Build up the indices
        indices = [i for i in 1:n]
        #All combinations of them
        x = collect(combinations(indices))
        #Now we want to normalize them
        #Meaning we set the first value of each pair to be 0
        x_c = []
        for a in x
            push!(x_c, a.-a[1])
        end
        #Now we collect how many unique correlations there are as well as degeneracy
        corrs, degen = unique_elements_and_frequencies(x_c)
        P = 1
        for a in 1:length(corrs)
            P+= degen[a]*abs(sigma_general(corrs[a],t))
         end
        return P/2^n
    end
    
    function PP_corr(n,l,t)
        indices = vcat([i for i in 1:n], [i for i in 1+l:n+l])
        x = collect(combinations(indices))
        x_c = []
        for a in x
            push!(x_c, a.-a[1])
        end
        #Now we collect how many unique correlations there are as well as degeneracy
        corrs, degen = unique_elements_and_frequencies(x_c)
        P = 1
        for a in 1:length(corrs)
            P+= degen[a]*abs(sigma_general(corrs[a],t))
        end
        return P/4^n- P_n(n,t)^2
    end
    
    
##Simple variables for testing only
J = parse(Float64,get(ENV, "J", "1"))
T_final = parse(Float64,get(ENV,"T_final","5"))
N_timesteps = parse(Int,get(ENV,"N_timesteps","10"))
L = parse(Int,get(ENV, "L", "50"))
h1 = parse(Float64,get(ENV, "h1","0.00"))
h2 = parse(Float64,get(ENV, "h2","1.50"))

println("J = $(J), L = $(L), h1 = $(h1), h2 = $(h2), N_t = $(N_timesteps)")
dt = (T_final-0.5)/N_timesteps
t = [i for i in 0.5:dt:T_final];
k = [2*pi*(n+1/2)/L for n in 0:L-1];
indices= []
plot_indices =[li for li in vcat(-10:1:-1,0,1:1:10)]
for n in 1:5
    l = [li for li in vcat(-10-n+1:1:-n,0,n:1:10+n-1)];
    push!(indices,l);
end
for ti in t
    PP_n = []
    for n in 1:5
        timer = time()
        push!(PP_n, [PP_corr(n,li,ti) for li in indices[n]])
        println("for n= $(n), time to compute = $(time()-timer)")
    end
    y = plot(plot_indices,PP_n[1],label="n=1")
    for n in 2:5
        plot!(plot_indices,PP_n[n],label="n= $(n)")
    end
    plot!(xlabel = "l")
    plot!(ylabel = "<P_nP_n(l)>")
    plot!(title = "h1 = $(h1), h2 = $(h2), time = $(ti)")
    savefig(y, "Waterfall_h1_$(h1)_h2_$(h2)_t_$(ti).png")
#############
#    y = plot(plot_indices,log.(PP_n[1]),label="n=1")
#    for n in 2:5
#        plot!(plot_indices,log.(PP_n[n]),label="n= $(n)")
#    end
#    plot!(xlabel = "l")
#    plot!(ylabel = "<P_nP_n(l)>")
#    plot!(title = "h1 = $(h1), h2 = $(h2), time = $(ti)")
#    savefig(y, "LOG_Waterfall_h1_$(h1)_h2_$(h2)_t_$(ti).png")

end

println("Succesfuly Completed")