
function u_tilde(τ,t,q)
    βq = 
    αq = 
    θ = exp(-1im 3* π/4)

    prefactor = θ* βq * sqrt(τ) * exp(-π βq^2 *τ/4 )
    nu = 1im*βq^2*τ-1
    arg = 1/θ * 2* (t-αq * τ)/sqrt(τ)

    return
end


using FewSpecialFunctions
V(1im+0.,1.0+1im)