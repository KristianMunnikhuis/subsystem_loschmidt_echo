using LinearAlgebra
using DifferentialEquations
using Combinatorics
function G_tfim(psi)
    return psi * psi'
end

function h_t(t, tau, h0, hf)
    return h0 + (hf - h0) * (t / tau)
end

function AA(Gi)
    L = size(Gi,1) ÷ 2
    G = Gi[1:L,1:L]
    F = Gi[1:L,L+1:end]
    return G + (I - G') + F + F'
end

function BB(Gi)
    L = size(Gi,1) ÷ 2
    G = Gi[1:L,1:L]
    F = Gi[1:L,L+1:end]
    return -G - (I - G') + F + F'
end

function AB(Gi)
    L = size(Gi,1) ÷ 2
    G = Gi[1:L,1:L]
    F = Gi[1:L,L+1:end]
    return G - (I - G') - F + F'
end

function BA(Gi)
    L = size(Gi,1) ÷ 2
    G = Gi[1:L,1:L]
    F = Gi[1:L,L+1:end]
    return -G + (I - G') - F + F'
end

function remove_duplicates_in_pairs(vec)
    unique_vals = unique(vec)
    counts = [count(==(x), vec) for x in unique_vals]
    return unique_vals[counts .% 2 .!= 0]
end

function sigma_general(indices, Gi, L)
    if length(indices)==0
        return 1
    end
    F = Gi[1:L,L+1:end]
    G = Gi[1:L,1:L]
    M = I - 2*(G+F)

    indices = sort(indices)
    indices = remove_duplicates_in_pairs(indices)
    
    if length(indices) % 2 == 1
        constant = 10
        indices = vcat(indices, indices .+ constant)
        return sqrt(abs(sigma_general(indices, Gi, L)))
    end

    odd_sites = indices[1:2:end]
    even_sites = indices[2:2:end]
    JW_string_lengths = even_sites .- odd_sites
    print(JW_string_lengths)
    N = sum(JW_string_lengths)

    R = Int[]
    for i in 1:2:length(indices)
        append!(R, odd_sites[i÷2+1]:even_sites[i÷2+1])
    end

    A_coords = setdiff(R, odd_sites)
    B_coords = setdiff(R, even_sites)

    C = zeros(ComplexF64, N, N)
    for nx in 1:N
        for ny in 1:N
            Bx = B_coords[nx]
            Ay = A_coords[ny]
            C[nx, ny] = M[Bx+1, Ay+1]
        end
    end

    return det(C)
end

function all_combinations(indices)
    result = [[]]
    for r in 1:length(indices)
        result = vcat(result, collect(combinations(indices, r)))
    end
    return result
end

function P_n(n, G, L)
    indices = collect(0:n-1)
    terms = all_combinations(indices)
    dat = [sigma_general(term, G, L) for term in terms]
    return mean(dat)
end

function H_bdg(h, L, J; boundary_condition="ABC")
    A = zeros(L, L)
    B = zeros(L, L)

    for j in 1:L
        A[j,j] = 2*h
    end
    for j in 1:(L-1)
        A[j,j+1] = A[j+1,j] = -J
        B[j,j+1] = -J
        B[j+1,j] = J
    end

    if boundary_condition == "ABC"
        A[1,L] = A[L,1] = J
        B[L,1] = J
        B[1,L] = -J
    elseif boundary_condition == "PBC"
        A[1,L] = A[L,1] = -J
        B[L,1] = -J
        B[1,L] = J
    end

    return (1/2) * [A B; -B -A]
end

function TFIM_time_evolve(N_steps, tau, h0, hf, J, L; U0=nothing, bc="ABC")
    function rhs!(dU, U, p, t)
        h = h_t(t, tau, h0, hf)
        H = H_bdg(h, L, J, boundary_condition=bc)
        Umat = reshape(U, 2L, 2L)
        dUmat = -im*2*H*Umat
        dU[:] = vec(dUmat)
        return nothing
    end

    if U0 === nothing
        H0 = H_bdg(h0, L, J, boundary_condition=bc)
        _, U0vals = eigen(H0)
        U0 = U0vals
    end
    U0vec = vec(complex(U0))

    tspan = (0.0, tau)
    tsteps = range(0.0, tau, length=N_steps)
    prob = ODEProblem(rhs!, U0vec, tspan)
    sol = solve(prob, Tsit5(), saveat=tsteps)

    U_t = [reshape(sol.u[i], 2L, 2L) for i in 1:length(sol.u)]
    return U_t
end
