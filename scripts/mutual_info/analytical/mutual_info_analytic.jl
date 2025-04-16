@doc """Analytical calculations of TFIM quantities"""

###Data Reading and Writing
base_folder = dirname(@__DIR__)
include( joinpath(@__DIR__,"..","..","..","src/analytical_TFIM.jl"))
plot_folder = joinpath(@__DIR__,"..","..","..","results/mutual_info/analytical/")
#mkdir(plot_folder)
#TFIM Functions 

x=  [sigma_general([1,2,3,4],ti) for ti in t]
y = sqrt.(abs.([sigma_general([1],ti) for ti in t]))

y[10:50] *= -1
plot(y)
plot!(sqrt.(real.(x)))
####PREFACTORS NEED WORK!!!!!
P_n(4,0)
P_n(6,0)

#Parameters
L = 40
J = 1
h1 = 0
h2 = 0
h_i = [i for i in 0:.01:2]
mut_1 = []
mut_2 = []
mut_3 = []
mut_4 = []
for hi in h_i[1:end]
    h1 = hi
    push!(mut_1,log(P_n(2,0.0)/P_n(1,0.0)^2))
    push!(mut_2,log(P_n(4,0.0)/P_n(2,0.0)^2))
    push!(mut_3,log(P_n(6,0.0)/P_n(3,0.0)^2))
    push!(mut_4,log(P_n(8,0.0)/P_n(4,0.0)^2))
end

plot(h_i,abs.(mut_1),dpi=300)
plot!(h_i,abs.(mut_2))
plot!(h_i,abs.(mut_3))
plot!(h_i,abs.(mut_4))
#Used for plotting
T_final = 5
t = [i for i in 0:0.05:T_final];
k = [2*pi*(n+1/2)/L for n in 0:L-1];

timer = time()
P10 = [[log(P_n(2*ni,ti)/(P_n(ni,ti)^2)) for ti in t] for ni in 1:3];

print(timer-time())
plot(t,real.(P10))

